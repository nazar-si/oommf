/* FILE: DMI_Cnv_z.cc            -*-Mode: c++-*-
 *
 * Dzyaloshinskii-Moriya energy for the Cnv crystallographic class [1]:
 *
 * $w_\text{dmi} = D ( L_{xz}^{(x)} + L_{yz}^{(y)} )
 *
 * This extension works both with and without periodic boundary conditions.
 *
 * Extension and modification by David Cortes-Ortuno, Marijan Beg and Hans Fangohr (University of Southampton and European XFEL GmbH) of Oxs_DMexchange6ngbr.h [2].
 *
 * Developed as a part of OpenDreamKit Horizon 2020 European Research Infrastructure
 * project (676541), and the EPSRC Programme grant on Skyrmionics (EP/N032128/1).
 *
 * [1] A. N. Bogdanov and D. A. Yablonskii. Zh. Eksp. Teor. Fiz. 95, 178-182 (1989).
 * [2] Rohart, S., & Thiaville, A. Physical Review B, 88, 184422 (2013).
 *
 */

#include <string>

#include "nb.h"
#include "key.h"
#include "director.h"
#include "mesh.h"
#include "meshvalue.h"
#include "simstate.h"
#include "threevector.h"
#include "rectangularmesh.h"
#include "DMI_Cnv_z.h"
#include "energy.h"		// Needed to make MSVC++ 5 happy

OC_USE_STRING;

// Oxs_Ext registration support
OXS_EXT_REGISTER(Oxs_DMI_Cnv_z);

/* End includes */


// Constructor
Oxs_DMI_Cnv_z::Oxs_DMI_Cnv_z(
  const char* name,     // Child instance id
  Oxs_Director* newdtr, // App director
  const char* argstr)   // MIF input block parameters
  : Oxs_Energy(name,newdtr,argstr),
    xperiodic(0), yperiodic(0), zperiodic(0),
    mesh_id(0)
{
  // Process arguments
  OXS_GET_INIT_EXT_OBJECT("D",Oxs_ScalarField, D_init)

  VerifyAllInitArgsUsed();
}

OC_BOOL Oxs_DMI_Cnv_z::Init()
{
  mesh_id = 0;
  return Oxs_Energy::Init();
}

void Oxs_DMI_Cnv_z::GetEnergy
(const Oxs_SimState& state,
 Oxs_EnergyData& oed
 ) const
{
  // See if mesh has changed.
  if(mesh_id !=  state.mesh->Id()) {
    // Setup region mapping
    mesh_id = 0; // Safety
    D_init->FillMeshValue(state.mesh,D);
    mesh_id = state.mesh->Id();
  }
  const Oxs_MeshValue<ThreeVector>& spin = state.spin;
  const Oxs_MeshValue<OC_REAL8m>& Ms_inverse = *(state.Ms_inverse);

  // Use supplied buffer space, and reflect that use in oed.
  oed.energy = oed.energy_buffer;
  oed.field = oed.field_buffer;
  Oxs_MeshValue<OC_REAL8m>& energy = *oed.energy_buffer;
  Oxs_MeshValue<ThreeVector>& field = *oed.field_buffer;

  // Check periodicity --------------------------------------------------------
  const Oxs_CommonRectangularMesh* mesh
    = dynamic_cast<const Oxs_CommonRectangularMesh*>(state.mesh);
  if(mesh==NULL) {
    String msg=String("Object ")
      + String(state.mesh->InstanceName())
      + String(" is not a rectangular mesh.");
    throw Oxs_ExtError(this,msg);
  }

  const Oxs_RectangularMesh* rmesh
    = dynamic_cast<const Oxs_RectangularMesh*>(mesh);
  const Oxs_PeriodicRectangularMesh* pmesh
    = dynamic_cast<const Oxs_PeriodicRectangularMesh*>(mesh);
  if(pmesh!=NULL) {
    // Rectangular, periodic mesh
    xperiodic = pmesh->IsPeriodicX();
    yperiodic = pmesh->IsPeriodicY();
    zperiodic = pmesh->IsPeriodicZ();
  } else if (rmesh!=NULL) {
    xperiodic=0; yperiodic=0; zperiodic=0;
  } else {
    String msg=String("Unknown mesh type: \"")
      + String(ClassName())
      + String("\".");
    throw Oxs_ExtError(this,msg.c_str());
  }
  // --------------------------------------------------------------------------

  OC_INDEX xdim = mesh->DimX();
  OC_INDEX ydim = mesh->DimY();
  OC_INDEX zdim = mesh->DimZ();
  OC_INDEX xydim = xdim * ydim;
  // OC_INDEX xyzdim = xdim * ydim * zdim;

  OC_REAL8m wgtx = 1.0/(mesh->EdgeLengthX());
  OC_REAL8m wgty = 1.0/(mesh->EdgeLengthY());
  //OC_REAL8m wgtz = -1.0/(mesh->EdgeLengthZ()*mesh->EdgeLengthZ());

  OC_REAL8m hcoef = -2/MU0;

  for(OC_INDEX z=0;z<zdim;z++) {
    for(OC_INDEX y=0;y<ydim;y++) {
      for(OC_INDEX x=0;x<xdim;x++) {
        OC_INDEX i = mesh->Index(x,y,z); // Get base linear address
        ThreeVector base = spin[i];
        OC_REAL8m Msii = Ms_inverse[i];
        if(Msii == 0.0) {
          energy[i]=0.0;
          field[i].Set(0.,0.,0.);
          continue;
        }
        ThreeVector sum(0.,0.,0.);
        ThreeVector zu(0.,0.,1.);
        OC_INDEX j;

        if(y > 0 || yperiodic) {  // y- direction
          if(y > 0) {
            j = i - xdim;
          } else if (yperiodic) {
            j = i - xdim + xydim;
          }
          if(Ms_inverse[j] != 0.0) {
            ThreeVector uij(0.,-1.,0);
            sum += 0.5 * D[i] * wgty * ((zu ^ uij) ^ spin[j]);
          }
        }

        if(x > 0 || xperiodic) {  // x- direction
          if(x > 0) {
            j = i - 1;        // j = mesh->Index(x-1,y,z)
          } else if (xperiodic) {
            j = i - 1 + xdim; // x == 0, j = Index(xdim-1,y,z);
          }
          if(Ms_inverse[j] != 0.0) {
            ThreeVector uij(-1.,0.,0);
            sum += 0.5 * D[i] * wgtx * ((zu ^ uij) ^ spin[j]);
          }
        }

        if(y < ydim - 1 || yperiodic) {  // y+ direction
          if (y < ydim-1) {
            j = i + xdim;
          } else if (yperiodic) {
            j = i + xdim - xydim;
          }
          if(Ms_inverse[j] != 0.0) {
            ThreeVector uij(0.,1.,0);
            sum += 0.5 * D[i] * wgty * ((zu ^ uij) ^ spin[j]);
          }
        }

        if(x < xdim-1 || xperiodic) {  // x+ direction
          if (x < xdim-1) {
            j = i + 1;
          } else if (xperiodic) {
            j = i + 1 - xdim;
          }
          if (Ms_inverse[j] != 0.0) {
            ThreeVector uij(1.,0.,0);
            sum += 0.5 * D[i] * wgtx * ((zu ^ uij) ^ spin[j]);
          }
        }
	
	    field[i]  = (hcoef * Msii) * sum;
	    energy[i] = (sum * base);
      }
    }
  }
}
