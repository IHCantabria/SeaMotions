
// Include local modules
#include "../math/integration.hpp"
#include "hmf_interface.hpp"


HMFInterface::HMFInterface(
                                SourceNode**    source_nodes,
                                cuscomplex*     source_values,
                                int             source_nodes_np,
                                PanelGeom*      panel,
                                int             start_index_i,
                                int             dof_j,
                                cusfloat        ang_freq,
                                cusfloat        water_depth,
                                cusfloat        grav_acc
                            )
{
    // Storage necessary input class arguments into the attributes
    this->_ang_freq         = ang_freq;
    this->_grav_acc         = grav_acc;
    this->_dof_j            = dof_j;
    this->_panel            = panel;
    this->_source_nodes     = source_nodes;
    this->_source_nodes_np  = source_nodes_np;
    this->_source_values    = source_values;
    this->_start_index_i    = start_index_i;
    this->_water_depth      = water_depth;

    // Generate green function interface
    this->_green_interf     = new   GWFInterface(
                                                    this->_source_nodes[0],
                                                    this->_source_values[0],
                                                    this->_panel->center,
                                                    this->_ang_freq,
                                                    this->_water_depth,
                                                    this->_grav_acc
                                                );

}


HMFInterface::~HMFInterface(
                                void
                            )
{
    delete this->_green_interf;
}


void        HMFInterface::set_ang_freq(
                                            cusfloat ang_freq
                                        )
{
    this->_ang_freq = ang_freq;
    this->_green_interf->set_ang_freq( ang_freq );
}


cuscomplex  HMFInterface::operator()(
                                        cusfloat ,
                                        cusfloat ,
                                        cusfloat x,
                                        cusfloat y,
                                        cusfloat z
                                    )
{
    // Create GaussPoint instatnce for the integration
    GaussPoints gp( 1 );

    // Set new field point for the integration
    cusfloat _field_point[3] = { x, y, z};
    this->_green_interf->set_field_point( _field_point );

    // Create lambda function for the integration of
    // the pressure over the panels
    GWFInterface*   gint    =   this->_green_interf;
    auto            lmb_fcn =   [gint]
                                (cusfloat xi, cusfloat eta, cusfloat x, cusfloat y, cusfloat z) -> cuscomplex
                                {
                                    return (*gint)( xi, eta, x, y, z );
                                };

    // Loop over source nodes list to integrate 
    // the pressure profile over the panel
    int         index       = 0;
    cuscomplex  potential   = std::complex( 0.0, 0.0 );
    cuscomplex  pressure    = std::complex( 0.0, 0.0 );
    for ( int i=0; i<this->_source_nodes_np; i++ )
    {
        // Set new source values
        index   =   this->_start_index_i + i;
        this->_green_interf->set_source(
                                            this->_source_nodes[i],
                                            this->_source_values[index]
                                        );

        // Integrate source value over the panel
        potential   +=  adaptive_quadrature_panel(
                                                    this->_panel,
                                                    lmb_fcn,
                                                    1000.0,
                                                    &gp
                                                );

        // Project pressure over the required direction
        pressure    +=  potential * this->_panel->normal_vec[this->_dof_j];
        
    }

    return pressure;
}


void    HMFInterface::set_start_index_i(
                                            int start_index
                                        )
{
    this->_start_index_i = start_index;
}


void    HMFInterface::set_dof_j(
                                    int dof_j
                                )
{
    this->_dof_j = dof_j;
}


void    HMFInterface::set_panel(
                                    PanelGeom* panel
                                )
{
    this->_panel = panel;
}


void    HMFInterface::set_source_values(
                                            cuscomplex* source_values
                                        )
{
    this->_source_values = source_values;
}