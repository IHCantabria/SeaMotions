
#ifndef __integration_math_hpp
#define __integration_math_hpp

// Include general usage libraries
#include <functional>

// Include local modules
#include "../config.hpp"
#include "../containers/panel_geom.hpp"
#include "../math/gauss.hpp"

// Generate alias for long type definitions in oder to clarify
// the function declaration
typedef std::function <cuscomplex(cuscomplex)> fcc_type;


/********************************************************/
/*************** Function Declaration *******************/
/********************************************************/
template<typename T>    inline  cuscomplex  adaptive_quadrature_panel(
                                                                            PanelGeom*      panel,
                                                                            T               fcn,
                                                                            cusfloat        tol,
                                                                            int             gp_order
                                                                        );

template<typename T>    inline  cuscomplex  adaptive_quadrature_panel(
                                                                            PanelGeom*      panel,
                                                                            T               fcn,
                                                                            cusfloat        tol,
                                                                            GaussPoints*    gp
                                                                        );

                                cuscomplex  complex_integration(
                                                                            fcc_type        f_def,
                                                                            cuscomplex      a,
                                                                            cuscomplex      b,
                                                                            cusfloat        tol
                                                                );    
                                            
                                cuscomplex  complex_integration_path(     
                                                                            fcc_type        f_def,
                                                                            int             num_segments,
                                                                            cuscomplex*     waypoints,
                                                                            cusfloat        tol,
                                                                            bool            close_path,
                                                                            bool            verbose
                                                                    );

template<typename T>            cuscomplex  quadrature_panel(
                                                                            PanelGeom*  panel,
                                                                            T           target_fcn,
                                                                            int         gp_order
                                                            );

template<typename T>            cuscomplex  quadrature_panel(
                                                                            PanelGeom*      panel,
                                                                            T               target_fcn,
                                                                            GaussPoints*    gp
                                                            );

template<typename Functor>      cusfloat    romberg_quadrature(   
                                                                            Functor         f, 
                                                                            cusfloat        a, 
                                                                            cusfloat        b, 
                                                                            cusfloat        precision
                                                                );

// Include templates definition
#include "integration.txx"

#endif