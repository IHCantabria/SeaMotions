
#ifndef __matlin_group_hpp
#define __matlin_group_hpp

// Include local modules
#include "../config.hpp"


template <typename T>
struct MatLinGroup
{
private:
    // Define class private attributes
    int         _dims_np            = 3;

public:
    // Define class public attributes
    cusfloat*   cog_to_field_points = nullptr;
    int         end_col             = 0;
    int         end_row             = 0;
    cusfloat*   field_points        = nullptr;
    int*        field_points_cnp    = nullptr;
    int         field_points_nb     = 0;
    int         field_points_np     = 0;
    T*          field_values        = nullptr;
    int         fields_np           = 0;
    int         start_col           = 0;
    int         start_row           = 0;
    T*          sysmat              = nullptr;
    int         sysmat_ncols        = 0;
    int         sysmat_nrows        = 0;
    T*          sysmat_steady       = nullptr;

    // Define class constructors and destructor
    MatLinGroup(
                    int sysmat_nrows_in,
                    int sysmat_ncols_in,
                    int field_points_nb_in,
                    int fields_np_in,
                    int start_row_in,
                    int end_row_in,
                    int start_col_in,
                    int end_col_in
                );

    ~MatLinGroup(
                    void
                );

    // Define class methods

};


template <typename T>
MatLinGroup<T>::MatLinGroup(
                                int sysmat_nrows_in,
                                int sysmat_ncols_in,
                                int field_points_nb_in,
                                int fields_np_in,
                                int start_row_in,
                                int end_row_in,
                                int start_col_in,
                                int end_col_in
                            )
{
    // Storage input arguments
    this->end_col           = end_col_in;
    this->end_row           = end_row_in;
    this->field_points_nb   = field_points_nb_in;
    this->field_points_np   = sysmat_nrows_in;
    this->fields_np         = fields_np_in;
    this->start_col         = start_col_in;
    this->start_row         = start_row_in;
    this->sysmat_nrows      = sysmat_nrows_in;
    this->sysmat_ncols      = sysmat_ncols_in;

    // Allocate space for the system matrixes
    this->sysmat            = generate_empty_vector<T>( this->sysmat_nrows * this->sysmat_ncols );
    this->sysmat_steady     = generate_empty_vector<T>( this->sysmat_nrows * this->sysmat_ncols );

    // Allocate space for the field points and the field values
    this->cog_to_field_points   = generate_empty_vector<cusfloat>( this->_dims_np * this->sysmat_nrows );
    this->field_points          = generate_empty_vector<cusfloat>( this->_dims_np * this->sysmat_nrows );
    this->field_points_cnp      = generate_empty_vector<int>( this->field_points_nb+1 );
    this->field_values          = generate_empty_vector<T>( this->fields_np * this->sysmat_nrows );
}


template <typename T>
MatLinGroup<T>::~MatLinGroup(
                            void
                        )
{
    mkl_free( this->cog_to_field_points );
    mkl_free( this->field_points );
    mkl_free( this->field_points_cnp );
    mkl_free( this->field_values );
    mkl_free( this->sysmat );
    mkl_free( this->sysmat_steady );
}


// Define short types to refer in a handly way to the MatLinGroup family members
typedef     MatLinGroup<cusfloat>       MLGFloat;
typedef     MatLinGroup<cuscomplex>     MLGCmpx;

#endif