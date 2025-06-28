
# Import general usage libraries
import os

# Import general usage libraries
import numpy as np

# Import general usage plotting libraries
from bokeh.plotting import figure, show, save, output_file
from bokeh.models import ColumnDataSource, LinearColorMapper, ColorBar, HoverTool
from bokeh.transform import transform
from bokeh.palettes import Viridis256


REF_COL = 60


class FitProperties:

    def __init__(self) -> None:
        self.cheby_order_x  = 0
        self.cheby_order_y  = 0
        self.cheby_order_z  = 0
        self.cheby_abs_tol  = 0.0
        self.fcn_log_scale  = False
        self.max_ref_level  = 0
        self.cheby_rel_tol  = 0.0
        self.num_x          = 60
        self.num_x_fit      = 30
        self.num_y          = 60
        self.num_y_fit      = 30
        self.num_z          = 60
        self.num_z_fit      = 30
        self.region_name    = ""
        self.x_log_scale    = False
        self.x_max          = 0
        self.x_map_fcn      = lambda x: x
        self.x_min          = 0.0
        self.y_log_scale    = False
        self.y_max          = 0.0
        self.y_map_fcn      = lambda y: y
        self.y_min          = 0.0
        self.z_log_scale    = False
        self.z_max          = 0.0
        self.z_map_fcn      = lambda z: z
        self.z_min          = 0.0

    def fit_points_to_order(self)->None:
        self.num_x = self.cheby_order_x
        self.num_x_fit = self.cheby_order_x
        self.num_y = self.cheby_order_y
        self.num_y_fit = self.cheby_order_y
        self.num_z = self.cheby_order_z
        self.num_z_fit = self.cheby_order_z

    def x_map(self, x: np.ndarray)->np.ndarray:
        return self.x_map_fcn(self.x_map_lin(x))
        # return self.x_map_fcn(x)

    def x_map_lin(self, x: np.ndarray)->np.ndarray:
        xmap = (x+1)*(self.x_max-self.x_min)/2 + self.x_min
        if self.x_log_scale:
            xmap = 10**xmap
        
        return xmap

    def y_map(self, y: np.ndarray)->np.ndarray:
        # return self.y_map_fcn(self.y_map_lin(y))
        return self.y_map_fcn(y)
    
    def y_map_lin(self, y: np.ndarray)->np.ndarray:
        ymap = (y+1)*(self.y_max-self.y_min)/2 + self.y_min
        if self.y_log_scale:
            ymap = 10**ymap
        
        return ymap

    def z_map(self, z: np.ndarray)->np.ndarray:
        # return self.y_map_fcn(self.y_map_lin(y))
        return self.z_map_fcn(z)
    
    def z_map_lin(self, z: np.ndarray)->np.ndarray:
        zmap = (z+1)*(self.z_max-self.z_min)/2 + self.z_min
        if self.z_log_scale:
            zmap = 10**zmap

        return zmap


class FitStats:

    def __init__(self) -> None:
        self.abs_err_over_thr   = 0.0
        self.abs_err_tol        = 0.0
        self.ref_level          = 0
        self.max_abs_err        = 0.0
        self.max_cheby_order    = 0
        self.max_rel_err        = 0.0
        self.max_value          = 0.0
        self.mean_abs_err       = 0.0
        self.mean_rel_err       = 0.0
        self.min_abs_err        = 0.0
        self.min_rel_err        = 0.0
        self.num_coeffs         = 0
        self.rel_err_over_thr   = 0.0
        self.rel_err_tol        = 0.0


class StatPatch:

    def __init__( self, name: str ) -> None:
        self.dx     = [ ]
        self.dy     = [ ]
        self.dz     = [ ]
        self.name   = name
        self.value  = [ ]
        self.x      = [ ]
        self.y      = [ ]
        self.z      = [ ]


class RefLevel:

    def __init__( self, fit_props: FitProperties, parent=None ) -> None:
        self.cheby_coeffs           = np.ndarray
        self.child                  = [ ]
        self.fit_props              = fit_props
        self.fit_stats              = FitProperties( )
        self.level                  = 0
        self.parent                 = parent
        self.start_index_global     = 0
        self.start_index_global_f   = 0

        if parent is not None:
            self.level = parent.level + 1

    def add_data( self, coeffs: np.ndarray, fit_stats: FitStats ) -> None:
        # Storage data
        self.cheby_coeffs           = coeffs
        self.fit_stats              = fit_stats

        # Set level value to fit_stats
        self.fit_stats.ref_level    = self.level

    def calculate_num_coeffs_folded( self ) -> int:
        if self.fit_props.dims == 1:
            ncf = 1

        elif self.fit_props.dims == 2:
            ncf = (
                        ( np.abs( np.diff( self.cheby_coeffs[:, 1] ) ) > 0.1 )
                    ).sum( ) + 1
            
        elif self.fit_props.dims == 3:
            ncf = (
                        ( np.abs( np.diff( self.cheby_coeffs[:, 1] ) ) > 0.1 )
                        +
                        ( np.abs( np.diff( self.cheby_coeffs[:, 2] ) ) > 0.1 )
                    ).sum( ) + 1

        return ncf

    def check_position_interval( self, pos: np.ndarray ):
        if pos.shape[0] == 1:
            is_pass =   ( 
                            ( self.fit_props.x_max > pos[0] ) & ( self.fit_props.x_min < pos[0] ) 
                        )
        elif pos.shape[0] == 2:
            is_pass =   ( 
                            ( self.fit_props.x_max > pos[0] ) & ( self.fit_props.x_min < pos[0] )
                            &
                            ( self.fit_props.y_max > pos[1] ) & ( self.fit_props.y_min < pos[1] )
                        )
            
        elif pos.shape[0] == 3:
            is_pass =   ( 
                            ( self.fit_props.x_max > pos[0] ) & ( self.fit_props.x_min < pos[0] )
                            &
                            ( self.fit_props.y_max > pos[1] ) & ( self.fit_props.y_min < pos[1] )
                            &
                            ( self.fit_props.z_max > pos[2] ) & ( self.fit_props.z_min < pos[2] )
                        )
        else:
            raise ValueError( "Could not check position in more that 3 dimensions." )
        
        return is_pass

    def check_tolerances( self ) -> bool:
        is_pass_abs         = self.fit_stats.max_abs_err < self.fit_stats.abs_err_tol
        is_pass_rel         = self.fit_stats.max_rel_err < self.fit_stats.rel_err_tol
        # is_pass_thr_abs     = self.fit_stats.abs_err_over_thr < 1.0
        is_pass_mean_abs    = self.fit_stats.mean_abs_err < self.fit_stats.abs_err_tol
        is_pass_loose_abs   = self.fit_stats.max_abs_err < 10 * self.fit_stats.abs_err_tol
        is_pass             = is_pass_abs or is_pass_rel or ( is_pass_mean_abs and is_pass_loose_abs )

        return is_pass
    
    def get_cheby_coeffs( self ) -> np.ndarray:
        if self.child:
            coeffs_i = [ ]
            for c in self.child:
                coeffs_i.append( c.get_cheby_coeffs( ) )

            coeffs_i = np.vstack( coeffs_i )

        else:
            coeffs_i = self.cheby_coeffs.copy( )

        return coeffs_i
    
    def calculate_max_cheby_order( self ) -> int:
        return self.cheby_coeffs[:, 1:].max( )
    
    def calculate_max_cheby_order_folded( self ) -> int:
        mcof = 0
        if self.fit_props.dims == 1:
            mcof = self.cheby_coeffs[:, 1].max( )
        
        elif self.fit_props.dims == 2:
            mcof = self.cheby_coeffs[:, 1:2].max( )

        elif self.fit_props.dims == 3:
            mcof = self.cheby_coeffs[:, 1:3].max( )

        return mcof
    
    def get_max_level( self ) -> int:
        if self.child:
            ref_level = -1
            for c in self.child:
                li = c.get_max_level( )
                if li > ref_level:
                    ref_level = li
        else:
            ref_level = self.level

        return ref_level

    def get_num_cheby_coeffs( self ) -> np.ndarray:
        if self.child:
            cum_coeffs  = [ ]
            last        = 0
            for c in self.child:
                for cc in c.get_num_cheby_coeffs( ):
                    last += cc
                    cum_coeffs.append( cc )

        else:
            cum_coeffs = [ self.cheby_coeffs.shape[0] ]

        return cum_coeffs
    
    def get_num_cheby_coeffs_folded( self ) -> np.ndarray:
        if self.child:
            cum_coeffs  = [ ]
            last        = 0
            for c in self.child:
                for cc in c.get_num_cheby_coeffs_folded( ):
                    last += cc
                    cum_coeffs.append( cc )

        else:
            cum_coeffs = [ self.calculate_num_coeffs_folded( ) ]

        return cum_coeffs

    def get_start_index( self, pos: np.ndarray ) -> int:
        if self.child:
            for c in self.child:
                index = c.get_start_index( pos )
                if index is not None:
                    break

        else:
            is_pos_interval = self.check_position_interval( pos )
            if is_pos_interval:
                index = (
                            self.start_index_global,
                            self.cheby_coeffs.shape[0],
                            self.start_index_global_f,
                            self.calculate_num_coeffs_folded( ),
                            self.calculate_max_cheby_order( ),
                            self.calculate_max_cheby_order_folded( ),
                            self.fit_props.x_min,
                            self.fit_props.x_max,
                            self.fit_props.y_min,
                            self.fit_props.y_max,
                            self.fit_props.z_min,
                            self.fit_props.z_max
                        )
                print( "Region -> ", index )
            
            else:
                index = None

        return index
    
    def get_stat_patch( self, stat_patch: StatPatch ) -> None:
        if self.child:
            for c in self.child:
                c.get_stat_patch( stat_patch )

        else:
            stat_patch.x.append( ( self.fit_props.x_max + self.fit_props.x_min ) / 2.0 )
            stat_patch.y.append( ( self.fit_props.y_max + self.fit_props.y_min ) / 2.0 )
            stat_patch.z.append( ( self.fit_props.z_max + self.fit_props.z_min ) / 2.0 )
            stat_patch.dx.append( self.fit_props.x_max - self.fit_props.x_min )
            stat_patch.dy.append( self.fit_props.y_max - self.fit_props.y_min )
            stat_patch.dz.append( self.fit_props.z_max - self.fit_props.z_min )
            stat_patch.value.append( getattr( self.fit_stats, stat_patch.name ) )

    def set_start_index( self, start_index: int ) -> int:
        if self.child:
            for c in self.child:
                start_index = c.set_start_index( start_index )
        
        else:
            self.start_index_global     = start_index + 0
            start_index                 += self.cheby_coeffs.shape[0]
            print( "Ref.Level: ", self.level, " - x_min: ", self.fit_props.x_min, " - x_max: ", self.fit_props.x_max, " - y_min: ", self.fit_props.y_min, " - y_max: ", self.fit_props.y_max, " - start_index: ", self.start_index_global )

        return start_index
    
    def set_start_index_folded( self, start_index: int ) -> int:
        if self.child:
            for c in self.child:
                start_index = c.set_start_index_folded( start_index )
        
        else:
            self.start_index_global_f   = start_index + 0
            start_index                 += self.calculate_num_coeffs_folded( )
            print( "Ref.Level: ", self.level, " - x_min: ", self.fit_props.x_min, " - x_max: ", self.fit_props.x_max, " - y_min: ", self.fit_props.y_min, " - y_max: ", self.fit_props.y_max, " - start_index: ", self.start_index_global )

        return start_index
    
    def show_summary( self, folder_path: str ) -> None:
        # Define output coefficients
        num_coeffs_sp       = StatPatch( "num_coeffs" )
        max_err_sp          = StatPatch( "max_abs_err" )
        ref_level_sp        = StatPatch( "ref_level" )
        max_cheby_order_sp  = StatPatch( "max_cheby_order" )
        max_abs_err_sp      = StatPatch( "max_abs_err" )
        mean_abs_err_sp     = StatPatch( "mean_abs_err" )
        min_abs_err_sp      = StatPatch( "min_abs_err" )

        self.get_stat_patch( num_coeffs_sp )
        self.get_stat_patch( max_err_sp )
        self.get_stat_patch( ref_level_sp )
        self.get_stat_patch( max_cheby_order_sp )
        self.get_stat_patch( max_abs_err_sp )
        self.get_stat_patch( mean_abs_err_sp )
        self.get_stat_patch( min_abs_err_sp )

        # Ouput summary data
        print( "Num.Coeffs: ", np.array( num_coeffs_sp.value ).sum( ) )
        print( "Num.Patches: ", len( max_err_sp.value ) )
        print( "Max.Cheyby Order: ", np.array( max_cheby_order_sp.value ).max( ) )
        with open( os.path.join( folder_path, "summary_stats.dat" ), "w" ) as fid:
            fid.writelines( f"Num.Coeffs:            {np.array( num_coeffs_sp.value ).sum( )}\n"  )
            fid.writelines( f"Num.Patches:           {len( max_err_sp.value )}\n" )
            fid.writelines( f"Max.Cheyby Order:      {np.array( max_cheby_order_sp.value ).max( )}" )

        plot_patch( num_coeffs_sp, folder_path )
        plot_patch( max_err_sp, folder_path )
        plot_patch( ref_level_sp, folder_path )
        plot_patch( max_cheby_order_sp, folder_path )
        plot_patch( max_abs_err_sp, folder_path )
        plot_patch( mean_abs_err_sp, folder_path )
        plot_patch( min_abs_err_sp, folder_path )


def plot_patch( stat_patch: StatPatch, folder_path: str ) -> None:
    # Create output file
    fipath = os.path.join( folder_path, stat_patch.name + ".html" )
    output_file( fipath, stat_patch.name.upper( ) )

    # Sample data
    x       = stat_patch.x
    y       = stat_patch.y
    width   = stat_patch.dx
    height  = stat_patch.dy
    value   = stat_patch.value

    # Prepare main rectangle source
    rect_data   = dict(x=x, y=y, width=width, height=height, value=value)
    rect_source = ColumnDataSource(rect_data)

    # Compute corner coordinates for each rectangle
    corner_xs = []
    corner_ys = []

    for xi, yi, wi, hi in zip(x, y, width, height):
        half_w = wi / 2
        half_h = hi / 2
        corners = [
            (xi - half_w, yi - half_h),
            (xi - half_w, yi + half_h),
            (xi + half_w, yi - half_h),
            (xi + half_w, yi + half_h),
        ]
        for cx, cy in corners:
            corner_xs.append(cx)
            corner_ys.append(cy)

    corner_source = ColumnDataSource(data=dict(x=corner_xs, y=corner_ys))

    # Set up color mapper
    color_mapper = LinearColorMapper(palette=Viridis256, low=min(value), high=max(value))

    # Create plot
    p = figure(
                    title=stat_patch.name.upper( ),
                    x_axis_label="X",
                    y_axis_label="Y",
                    tools="pan,box_zoom,wheel_zoom,reset",
                    sizing_mode="scale_height"
                )

    # Draw rectangles with hover
    rects = p.rect(x='x', y='y', width='width', height='height', 
                fill_color=transform('value', color_mapper), line_color='black', source=rect_source)

    # Add hover tool for rectangles
    hover = HoverTool(renderers=[rects],
                    tooltips=[
                        ("Value", "@value"),
                        ("Center (x, y)", "(@x, @y)")
                    ])
    p.add_tools(hover)

    # Draw corner circles
    p.circle(x='x', y='y', size=6, color='black', source=corner_source)

    # Add color bar
    color_bar = ColorBar(color_mapper=color_mapper, label_standoff=12, location=(0, 0))
    p.add_layout(color_bar, 'right')

    save( p )
    show( p )


def str_start_col(str0: str, str1: str, num_col: int)->str:
    l0 = len(str0)
    if l0 >= num_col:
        raise ValueError("First string is bigger than the reference column.")

    for i in range(num_col-l0):
        str0 += " "
    
    return str0+str1


def write_coeffs_module_adaptive_1d_only_header( ref_level: RefLevel, folder_path: str, module_name: str)->None:
    # Get number of cumulative coefficients
    num_points_cum = ref_level.get_num_cheby_coeffs( )

    # Get cumulative data
    cum_coeffs      = ref_level.get_cheby_coeffs( )
    cheby_coeffs    = cum_coeffs[:, 0]
    ncx             = cum_coeffs[:, 1]
    max_cheby_order = cum_coeffs[:, 1:].max( )

    # Get intervals size
    max_ref_level   = ref_level.get_max_level( )
    x_min           = ref_level.fit_props.x_min
    x_max           = ref_level.fit_props.x_max

    # Get log scales
    x_log_scale     = ref_level.fit_props.x_log_scale
    fcn_log_scale   = ref_level.fit_props.fcn_log_scale

    # Calculate hash table
    intervals_np    = 2**(max_ref_level)
    dx              = ( x_max - x_min ) / intervals_np
    x               = np.arange( x_min, x_max+dx, dx )
    xm              = ( x[1:] + x[:-1] ) / 2.0

    blocks_np                   = xm.shape[0]
    blocks_start                = np.zeros( ( blocks_np, ), dtype=int )
    blocks_coeffs_np            = np.zeros( ( blocks_np, ), dtype=int )
    blocks_start_f              = np.zeros( ( blocks_np, ), dtype=int )
    blocks_coeffs_np_f          = np.zeros( ( blocks_np, ), dtype=int )
    blocks_max_cheby_order      = np.zeros( ( blocks_np, ), dtype=int )
    blocks_max_cheby_order_f    = np.zeros( ( blocks_np, ), dtype=int )
    x_min_vec                   = np.zeros( ( blocks_np, ), dtype=float )
    x_max_vec                   = np.zeros( ( blocks_np, ), dtype=float )
    count                       = 0
    for i in range( xm.shape[0] ):
        xmi = xm[i]
        bs, bc, bsf, bcf, mco, mcof, x_min_i, x_max_i, y_min_i, y_max_i,_, _    = ref_level.get_start_index( np.array( [ xmi ] ) )
        blocks_start[count]                                                     = bs
        blocks_coeffs_np[count]                                                 = bc
        blocks_start_f[count]                                                   = bsf
        blocks_coeffs_np_f[count]                                               = bcf
        blocks_max_cheby_order[count]                                           = mco
        blocks_max_cheby_order_f[count]                                         = mcof
        x_min_vec[count]                                                        = x_min_i
        x_max_vec[count]                                                        = x_max_i
        count                                                                   += 1

    max_cheby_order_f = blocks_max_cheby_order_f.max( )
    dx_vec = x_max_vec - x_min_vec
    
    # Open file unit to storage the coefficients in the 
    # database
    header_file_path = os.path.join(folder_path, module_name+".hpp")
    fid = open(header_file_path, "w")

    # Add header guard
    fid.writelines(f"#ifndef __{module_name}_coeffs_hpp\n")
    fid.writelines(f"#define __{module_name}_coeffs_hpp\n\n")

    # Add includes
    fid.writelines('#include "../../config.hpp"\n')
    fid.writelines("\n")

    # Open namespace field
    fid.writelines(f"struct {module_name}" + "C\n{\n")

    # Save number of intervals
    fid.writelines(str_start_col(f"static constexpr int", f"max_ref_level = {max_ref_level:d};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr int", f"intervals_np = {intervals_np:d};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr int", f"max_cheby_order = {int(max_cheby_order):d};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr int", f"max_cheby_order_f = {int(max_cheby_order_f):d};\n", REF_COL))

    # Save number of points
    blocks_start_str = ", ".join(f"{i:d}" for i in blocks_start)
    fid.writelines(str_start_col(f"static constexpr std::size_t", f"blocks_start[{blocks_np}] = " + "{" + blocks_start_str + "};\n", REF_COL))
    blocks_start_f_str = ", ".join(f"{i:d}" for i in blocks_start_f)
    fid.writelines(str_start_col(f"static constexpr std::size_t", f"blocks_start_fall[{blocks_np}] = " + "{" + blocks_start_f_str + "};\n", REF_COL))
    blocks_coeffs_np_str = ", ".join(f"{i:d}" for i in blocks_coeffs_np)
    fid.writelines(str_start_col(f"static constexpr std::size_t", f"blocks_coeffs_np[{blocks_np}] = " + "{" + blocks_coeffs_np_str + "};\n", REF_COL))
    blocks_coeffs_np_f_str = ", ".join(f"{i:d}" for i in blocks_coeffs_np_f)
    fid.writelines(str_start_col(f"static constexpr std::size_t", f"blocks_coeffs_np_fall[{blocks_np}] = " + "{" + blocks_coeffs_np_f_str + "};\n", REF_COL))
    blocks_max_cheby_order_str = ", ".join(f"{i:d}" for i in blocks_max_cheby_order)
    fid.writelines(str_start_col(f"static constexpr std::size_t", f"blocks_max_cheby_order[{blocks_np}] = " + "{" + blocks_max_cheby_order_str + "};\n", REF_COL))
    blocks_max_cheby_order_f_str = ", ".join(f"{i:d}" for i in blocks_max_cheby_order_f)
    fid.writelines(str_start_col(f"static constexpr std::size_t", f"blocks_max_cheby_order_fall[{blocks_np}] = " + "{" + blocks_max_cheby_order_f_str + "};\n", REF_COL))
    x_min_str = ", ".join(f"{v:0.3E}" for v in x_min_vec)
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"x_min_region[{blocks_np}] = " + "{" + x_min_str + "};\n", REF_COL))
    x_max_str = ", ".join(f"{v:0.3E}" for v in x_max_vec)
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"x_max_region[{blocks_np}] = " + "{" + x_max_str + "};\n", REF_COL))
    dx_str = ", ".join(f"{v:0.3E}" for v in dx_vec)
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"dx_region[{blocks_np}] = " + "{" + dx_str + "};\n", REF_COL))

    # Save interval bounds
    fid.writelines(str_start_col(f"static constexpr bool", f"fcn_log_scale = {fcn_log_scale:d};\n", REF_COL))

    fid.writelines(str_start_col(f"static constexpr bool", f"x_log_scale = {x_log_scale:d};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"x_max_global = {x_max:0.16f};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"x_min_global = {x_min:0.16f};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"dx_min_region = {dx:0.16f};\n", REF_COL))

    fid.writelines("\n")

    # Write chebyshev polynomials
    fid.writelines(str_start_col("static constexpr std::size_t", f"num_c = {cheby_coeffs.shape[0]};\n", REF_COL))
    write_vector_header(fid, cheby_coeffs, "c", "cusfloat")

    # Write polynomials coefficients
    write_vector_header(fid, ncx, "ncx", "std::size_t")

    # Close namespace field
    fid.writelines("};\n")

    # Close header guard if statement
    fid.writelines("#endif\n")

    # Close file unit for header file
    fid.close()


def write_coeffs_module_adaptive_2d_only_header( ref_level: RefLevel, folder_path: str, module_name: str)->None:
    # Get number of cumulative coefficients
    num_points_cum = ref_level.get_num_cheby_coeffs( )

    # Get cumulative data
    cum_coeffs      = ref_level.get_cheby_coeffs( )
    cheby_coeffs    = cum_coeffs[:, 0]
    ncx             = cum_coeffs[:, 1]
    ncy             = cum_coeffs[:, 2]
    max_cheby_order = cum_coeffs[:, 1:].max( )

    # Get intervals size
    max_ref_level   = ref_level.get_max_level( )
    x_min           = ref_level.fit_props.x_min
    x_max           = ref_level.fit_props.x_max
    y_min           = ref_level.fit_props.y_min
    y_max           = ref_level.fit_props.y_max

    # Get log scales
    x_log_scale     = ref_level.fit_props.x_log_scale
    y_log_scale     = ref_level.fit_props.y_log_scale
    fcn_log_scale   = ref_level.fit_props.fcn_log_scale

    # Calculate hash table
    intervals_np    = 2**(max_ref_level)
    dx              = ( x_max - x_min ) / intervals_np
    x               = np.arange( x_min, x_max+dx, dx )
    xm              = ( x[1:] + x[:-1] ) / 2.0

    dy              = ( y_max - y_min ) / intervals_np
    y               = np.arange( y_min, y_max+dy, dy )
    ym              = ( y[1:] + y[:-1] ) / 2.0

    blocks_np                   = xm.shape[0] * ym.shape[0]
    blocks_start                = np.zeros( ( blocks_np, ), dtype=int )
    blocks_coeffs_np            = np.zeros( ( blocks_np, ), dtype=int )
    blocks_max_cheby_order      = np.zeros( ( blocks_np, ), dtype=int )
    x_min_vec                   = np.zeros( ( blocks_np, ), dtype=float )
    x_max_vec                   = np.zeros( ( blocks_np, ), dtype=float )
    y_min_vec                   = np.zeros( ( blocks_np, ), dtype=float )
    y_max_vec                   = np.zeros( ( blocks_np, ), dtype=float )
    count                       = 0
    for i in range( xm.shape[0] ):
        xmi = xm[i]
        for j in range( ym.shape[0] ):
            ymi = ym[j]
            bs, bc, _, _, mco, _, x_min_i, x_max_i, y_min_i, y_max_i, _, _      = ref_level.get_start_index( np.array( [ xmi, ymi ] ) )
            blocks_start[count]                                                 = bs
            blocks_coeffs_np[count]                                             = bc
            blocks_max_cheby_order[count]                                       = mco
            x_min_vec[count]                                                    = x_min_i
            x_max_vec[count]                                                    = x_max_i
            y_min_vec[count]                                                    = y_min_i
            y_max_vec[count]                                                    = y_max_i
            count                                                               += 1

    dx_vec = x_max_vec - x_min_vec
    dy_vec = y_max_vec - y_min_vec
        
    # Open file unit to storage the coefficients in the 
    # database
    header_file_path = os.path.join(folder_path, module_name+".hpp")
    fid = open(header_file_path, "w")

    # Add header guard
    fid.writelines(f"#ifndef __{module_name}_coeffs_hpp\n")
    fid.writelines(f"#define __{module_name}_coeffs_hpp\n\n")

    # Add includes
    fid.writelines('#include "../../config.hpp"\n')
    fid.writelines("\n")

    # Open namespace field
    fid.writelines(f"struct {module_name}" + "C\n{\n")

    # Save number of intervals
    fid.writelines(str_start_col(f"static constexpr int", f"max_ref_level = {max_ref_level:d};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr int", f"intervals_np = {intervals_np:d};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr int", f"max_cheby_order = {int(max_cheby_order):d};\n", REF_COL))

    # Save number of points
    blocks_start_str = ", ".join(f"{i:d}" for i in blocks_start)
    fid.writelines(str_start_col(f"static constexpr std::size_t", f"blocks_start[{blocks_np}] = " + "{" + blocks_start_str + "};\n", REF_COL))
    blocks_coeffs_np_str = ", ".join(f"{i:d}" for i in blocks_coeffs_np)
    fid.writelines(str_start_col(f"static constexpr std::size_t", f"blocks_coeffs_np[{blocks_np}] = " + "{" + blocks_coeffs_np_str + "};\n", REF_COL))
    blocks_max_cheby_order_str = ", ".join(f"{i:d}" for i in blocks_max_cheby_order)
    fid.writelines(str_start_col(f"static constexpr std::size_t", f"blocks_max_cheby_order[{blocks_np}] = " + "{" + blocks_max_cheby_order_str + "};\n", REF_COL))
    x_min_str = ", ".join(f"{v:0.3E}" for v in x_min_vec)
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"x_min_region[{blocks_np}] = " + "{" + x_min_str + "};\n", REF_COL))
    x_max_str = ", ".join(f"{v:0.3E}" for v in x_max_vec)
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"x_max_region[{blocks_np}] = " + "{" + x_max_str + "};\n", REF_COL))
    dx_str = ", ".join(f"{v:0.3E}" for v in dx_vec)
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"dx_region[{blocks_np}] = " + "{" + dx_str + "};\n", REF_COL))
    y_min_str = ", ".join(f"{v:0.3E}" for v in y_min_vec)
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"y_min_region[{blocks_np}] = " + "{" + y_min_str + "};\n", REF_COL))
    y_max_str = ", ".join(f"{v:0.3E}" for v in y_max_vec)
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"y_max_region[{blocks_np}] = " + "{" + y_max_str + "};\n", REF_COL))
    dy_str = ", ".join(f"{v:0.3E}" for v in dy_vec)
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"dy_region[{blocks_np}] = " + "{" + dy_str + "};\n", REF_COL))

    # Save interval bounds
    fid.writelines(str_start_col(f"static constexpr bool", f"fcn_log_scale = {fcn_log_scale:d};\n", REF_COL))

    fid.writelines(str_start_col(f"static constexpr bool", f"x_log_scale = {x_log_scale:d};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"x_max_global = {x_max:0.16f};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"x_min_global = {x_min:0.16f};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"dx_min_region = {dx:0.16f};\n", REF_COL))

    fid.writelines(str_start_col(f"static constexpr bool", f"y_log_scale = {y_log_scale:d};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"y_max_global = {y_max:0.16f};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"y_min_global = {y_min:0.16f};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"dy_min_region = {dy:0.16f};\n", REF_COL))

    fid.writelines("\n")

    # Write chebyshev polynomials
    fid.writelines(str_start_col("static constexpr std::size_t", f"num_c = {cheby_coeffs.shape[0]};\n", REF_COL))
    write_vector_header(fid, cheby_coeffs, "c", "cusfloat")

    # Write polynomials coefficients
    write_vector_header(fid, ncx, "ncx", "std::size_t")
    write_vector_header(fid, ncy, "ncy", "std::size_t")


    # Close namespace field
    fid.writelines("};\n")

    # Close header guard if statement
    fid.writelines("#endif\n")

    # Close file unit for header file
    fid.close()


def write_coeffs_module_adaptive_3d_only_header( ref_level: RefLevel, folder_path: str, module_name: str)->None:
    # Get number of cumulative coefficients
    num_points_cum = ref_level.get_num_cheby_coeffs( )

    # Get cumulative data
    cum_coeffs      = ref_level.get_cheby_coeffs( )
    cheby_coeffs    = cum_coeffs[:, 0]
    ncx             = cum_coeffs[:, 1]
    ncy             = cum_coeffs[:, 2]
    ncz             = cum_coeffs[:, 3]
    ncf             = np.array( ref_level.get_num_cheby_coeffs_folded( ) ).sum( )
    max_cheby_order = cum_coeffs[:, 1:].max( )

    # Get intervals size
    max_ref_level   = ref_level.get_max_level( )
    x_min           = ref_level.fit_props.x_min
    x_max           = ref_level.fit_props.x_max
    y_min           = ref_level.fit_props.y_min
    y_max           = ref_level.fit_props.y_max
    z_min           = ref_level.fit_props.z_min
    z_max           = ref_level.fit_props.z_max

    # Get log scales
    x_log_scale     = ref_level.fit_props.x_log_scale
    y_log_scale     = ref_level.fit_props.y_log_scale
    z_log_scale     = ref_level.fit_props.z_log_scale
    fcn_log_scale   = ref_level.fit_props.fcn_log_scale

    # Calculate hash table
    intervals_np    = 2**(max_ref_level)
    dx              = ( x_max - x_min ) / intervals_np
    x               = np.arange( x_min, x_max+dx, dx )
    xm              = ( x[1:] + x[:-1] ) / 2.0

    dy              = ( y_max - y_min ) / intervals_np
    y               = np.arange( y_min, y_max+dy, dy )
    ym              = ( y[1:] + y[:-1] ) / 2.0

    dz              = ( z_max - z_min ) / intervals_np
    z               = np.arange( z_min, z_max+dz, dz )
    zm              = ( z[1:] + z[:-1] ) / 2.0

    blocks_np                   = xm.shape[0] * ym.shape[0] * zm.shape[0]
    blocks_np_f                 = xm.shape[0] * ym.shape[0]
    blocks_start                = np.zeros( ( blocks_np, ), dtype=int )
    blocks_start_f              = np.zeros( ( blocks_np, ), dtype=int )
    blocks_coeffs_np            = np.zeros( ( blocks_np, ), dtype=int )
    blocks_coeffs_np_f          = np.zeros( ( blocks_np, ), dtype=int )
    blocks_max_cheby_order      = np.zeros( ( blocks_np, ), dtype=int )
    blocks_max_cheby_order_f    = np.zeros( ( blocks_np, ), dtype=int )
    x_min_vec                   = np.zeros( ( blocks_np, ), dtype=float )
    x_max_vec                   = np.zeros( ( blocks_np, ), dtype=float )
    y_min_vec                   = np.zeros( ( blocks_np, ), dtype=float )
    y_max_vec                   = np.zeros( ( blocks_np, ), dtype=float )
    z_min_vec                   = np.zeros( ( blocks_np, ), dtype=float )
    z_max_vec                   = np.zeros( ( blocks_np, ), dtype=float )
    count                       = 0
    for i in range( xm.shape[0] ):
        xmi = xm[i]
        for j in range( ym.shape[0] ):
            ymi = ym[j]
            for k in range( zm.shape[0] ):
                zmi = zm[k]
                bs, bc, bsf, bcf, mco, mcof, x_min_i, x_max_i, y_min_i, y_max_i, z_min_i, z_max_i       = ref_level.get_start_index( np.array( [ xmi, ymi, zmi ] ) )
                blocks_start[count]                                                                     = bs
                blocks_coeffs_np[count]                                                                 = bc
                blocks_start_f[count]                                                                   = bsf
                blocks_coeffs_np_f[count]                                                               = bcf
                blocks_max_cheby_order[count]                                                           = mco
                blocks_max_cheby_order_f[count]                                                         = mcof
                x_min_vec[count]                                                                        = x_min_i
                x_max_vec[count]                                                                        = x_max_i
                y_min_vec[count]                                                                        = y_min_i
                y_max_vec[count]                                                                        = y_max_i
                z_min_vec[count]                                                                        = z_min_i
                z_max_vec[count]                                                                        = z_max_i
                count                                                                                   += 1

    max_cheby_order_f = blocks_max_cheby_order_f.max( )

    dx_vec = x_max_vec - x_min_vec
    dy_vec = y_max_vec - y_min_vec
    dz_vec = z_max_vec - z_min_vec
        
    # Open file unit to storage the coefficients in the 
    # database
    header_file_path = os.path.join(folder_path, module_name+".hpp")
    fid = open(header_file_path, "w")

    # Add header guard
    fid.writelines(f"#ifndef __{module_name}_coeffs_hpp\n")
    fid.writelines(f"#define __{module_name}_coeffs_hpp\n\n")

    # Add includes
    fid.writelines('#include "../../config.hpp"\n')
    fid.writelines("\n")

    # Open namespace field
    fid.writelines(f"struct {module_name}" + "C\n{\n")

    # Save number of intervals
    fid.writelines(str_start_col(f"static constexpr int", f"max_ref_level = {max_ref_level:d};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr int", f"intervals_np = {intervals_np:d};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr int", f"max_cheby_order = {int(max_cheby_order):d};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr int", f"max_cheby_order_f = {int(max_cheby_order_f):d};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr int", f"blocks_np = {int(blocks_np):d};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr int", f"blocks_np_f = {int(blocks_np_f):d};\n", REF_COL))

    # Save number of points
    blocks_start_str = ", ".join(f"{i:d}" for i in blocks_start)
    fid.writelines(str_start_col(f"static constexpr std::size_t", f"blocks_start[{blocks_np}] = " + "{" + blocks_start_str + "};\n", REF_COL))
    blocks_start_f_str = ", ".join(f"{i:d}" for i in blocks_start_f)
    fid.writelines(str_start_col(f"static constexpr std::size_t", f"blocks_start_fall[{blocks_np}] = " + "{" + blocks_start_f_str + "};\n", REF_COL))
    blocks_coeffs_np_str = ", ".join(f"{i:d}" for i in blocks_coeffs_np)
    fid.writelines(str_start_col(f"static constexpr std::size_t", f"blocks_coeffs_np[{blocks_np}] = " + "{" + blocks_coeffs_np_str + "};\n", REF_COL))
    blocks_coeffs_np_f_str = ", ".join(f"{i:d}" for i in blocks_coeffs_np_f)
    fid.writelines(str_start_col(f"static constexpr std::size_t", f"blocks_coeffs_np_fall[{blocks_np}] = " + "{" + blocks_coeffs_np_f_str + "};\n", REF_COL))
    blocks_max_cheby_order_str = ", ".join(f"{i:d}" for i in blocks_max_cheby_order)
    fid.writelines(str_start_col(f"static constexpr std::size_t", f"blocks_max_cheby_order[{blocks_np}] = " + "{" + blocks_max_cheby_order_str + "};\n", REF_COL))
    blocks_max_cheby_order_f_str = ", ".join(f"{i:d}" for i in blocks_max_cheby_order_f)
    fid.writelines(str_start_col(f"static constexpr std::size_t", f"blocks_max_cheby_order_fall[{blocks_np}] = " + "{" + blocks_max_cheby_order_f_str + "};\n", REF_COL))
    x_min_str = ", ".join(f"{v:0.3E}" for v in x_min_vec)
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"x_min_region[{blocks_np}] = " + "{" + x_min_str + "};\n", REF_COL))
    x_max_str = ", ".join(f"{v:0.3E}" for v in x_max_vec)
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"x_max_region[{blocks_np}] = " + "{" + x_max_str + "};\n", REF_COL))
    dx_str = ", ".join(f"{v:0.3E}" for v in dx_vec)
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"dx_region[{blocks_np}] = " + "{" + dx_str + "};\n", REF_COL))
    y_min_str = ", ".join(f"{v:0.3E}" for v in y_min_vec)
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"y_min_region[{blocks_np}] = " + "{" + y_min_str + "};\n", REF_COL))
    y_max_str = ", ".join(f"{v:0.3E}" for v in y_max_vec)
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"y_max_region[{blocks_np}] = " + "{" + y_max_str + "};\n", REF_COL))
    dy_str = ", ".join(f"{v:0.3E}" for v in dy_vec)
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"dy_region[{blocks_np}] = " + "{" + dy_str + "};\n", REF_COL))
    z_min_str = ", ".join(f"{v:0.3E}" for v in z_min_vec)
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"z_min_region[{blocks_np}] = " + "{" + z_min_str + "};\n", REF_COL))
    z_max_str = ", ".join(f"{v:0.3E}" for v in z_max_vec)
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"z_max_region[{blocks_np}] = " + "{" + z_max_str + "};\n", REF_COL))
    dz_str = ", ".join(f"{v:0.3E}" for v in dz_vec)
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"dz_region[{blocks_np}] = " + "{" + dz_str + "};\n", REF_COL))

    blocks_start_f_str = ", ".join(f"{0:d}" for _ in range( blocks_np_f ))
    fid.writelines(str_start_col(f"inline static std::size_t", f"blocks_start_f[{blocks_np_f}] = " + "{" + blocks_start_f_str + "};\n", REF_COL))
    blocks_coeffs_np_f_str = ", ".join(f"{0:d}" for _ in range( blocks_np_f ))
    fid.writelines(str_start_col(f"inline static std::size_t", f"blocks_coeffs_np_f[{blocks_np_f}] = " + "{" + blocks_coeffs_np_f_str + "};\n", REF_COL))
    blocks_max_cheby_order_f_str = ", ".join(f"{0:d}" for _ in range( blocks_np_f ))
    fid.writelines(str_start_col(f"inline static std::size_t", f"blocks_max_cheby_order_f[{blocks_np_f}] = " + "{" + blocks_max_cheby_order_f_str + "};\n", REF_COL))
    x_min_f_str = ", ".join(f"{0:0.3E}" for _ in range( blocks_np_f ))
    fid.writelines(str_start_col(f"inline static cusfloat", f"x_min_region_f[{blocks_np_f}] = " + "{" + x_min_f_str + "};\n", REF_COL))
    x_max_f_str = ", ".join(f"{0:0.3E}" for _ in range( blocks_np_f ))
    fid.writelines(str_start_col(f"inline static cusfloat", f"x_max_region_f[{blocks_np_f}] = " + "{" + x_max_f_str + "};\n", REF_COL))
    dx_f_str = ", ".join(f"{0:0.3E}" for _ in range( blocks_np_f ))
    fid.writelines(str_start_col(f"inline static cusfloat", f"dx_region_f[{blocks_np_f}] = " + "{" + dx_f_str + "};\n", REF_COL))
    y_min_f_str = ", ".join(f"{0:0.3E}" for _ in range( blocks_np_f ))
    fid.writelines(str_start_col(f"inline static cusfloat", f"y_min_region_f[{blocks_np_f}] = " + "{" + y_min_f_str + "};\n", REF_COL))
    y_max_f_str = ", ".join(f"{0:0.3E}" for _ in range( blocks_np_f ))
    fid.writelines(str_start_col(f"inline static cusfloat", f"y_max_region_f[{blocks_np_f}] = " + "{" + y_max_f_str + "};\n", REF_COL))
    dy_f_str = ", ".join(f"{0:0.3E}" for _ in range( blocks_np_f ))
    fid.writelines(str_start_col(f"inline static cusfloat", f"dy_region_f[{blocks_np_f}] = " + "{" + dy_f_str + "};\n", REF_COL))
    z_min_f_str = ", ".join(f"{0:0.3E}" for _ in range( blocks_np_f ))
    fid.writelines(str_start_col(f"inline static cusfloat", f"z_min_region_f[{blocks_np_f}] = " + "{" + z_min_f_str + "};\n", REF_COL))
    z_max_f_str = ", ".join(f"{0:0.3E}" for _ in range( blocks_np_f ))
    fid.writelines(str_start_col(f"inline static cusfloat", f"z_max_region_f[{blocks_np_f}] = " + "{" + z_max_f_str + "};\n", REF_COL))
    dz_f_str = ", ".join(f"{0:0.3E}" for _ in range( blocks_np_f ))
    fid.writelines(str_start_col(f"inline static cusfloat", f"dz_region_f[{blocks_np_f}] = " + "{" + dz_f_str + "};\n", REF_COL))

    # Save interval bounds
    fid.writelines(str_start_col(f"static constexpr bool", f"fcn_log_scale = {fcn_log_scale:d};\n", REF_COL))

    fid.writelines(str_start_col(f"static constexpr bool", f"x_log_scale = {x_log_scale:d};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"x_max_global = {x_max:0.16f};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"x_min_global = {x_min:0.16f};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"dx_min_region = {dx:0.16f};\n", REF_COL))

    fid.writelines(str_start_col(f"static constexpr bool", f"y_log_scale = {y_log_scale:d};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"y_max_global = {y_max:0.16f};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"y_min_global = {y_min:0.16f};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"dy_min_region = {dy:0.16f};\n", REF_COL))

    fid.writelines(str_start_col(f"static constexpr bool", f"z_log_scale = {z_log_scale:d};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"z_max_global = {z_max:0.16f};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"z_min_global = {z_min:0.16f};\n", REF_COL))
    fid.writelines(str_start_col(f"static constexpr cusfloat", f"dz_min_region = {dz:0.16f};\n", REF_COL))

    fid.writelines("\n")

    # Write chebyshev polynomials
    fid.writelines(str_start_col("static constexpr std::size_t", f"num_c = {cheby_coeffs.shape[0]};\n", REF_COL))
    fid.writelines(str_start_col("static constexpr std::size_t", f"num_cf = {ncf};\n", REF_COL))
    write_vector_header(fid, cheby_coeffs, "c", "cusfloat")
    write_vector_header(fid, np.zeros( ( ncf, ) ), "cf", "cusfloat", keyword="")

    # Write polynomials coefficients
    write_vector_header(fid, ncx, "ncx", "std::size_t")
    write_vector_header(fid, np.zeros( ( ncf, ), dtype=int ), "ncxf", "std::size_t", keyword="", is_inline=True)
    write_vector_header(fid, ncy, "ncy", "std::size_t")
    write_vector_header(fid, np.zeros( ( ncf, ), dtype=int ), "ncyf", "std::size_t", keyword="", is_inline=True)
    write_vector_header(fid, ncz, "ncz", "std::size_t")


    # Close namespace field
    fid.writelines("};\n")

    # Close header guard if statement
    fid.writelines("#endif\n")

    # Close file unit for header file
    fid.close()


def write_vector(fid, field: np.ndarray, field_tag: str, var_type: str, int_name: str)->None:
    fid.writelines(str_start_col(f"alignas(FLOATING_PRECISION) const {var_type}", f"{int_name}::{field_tag}[{field.shape[0]}] = " + "{\n", REF_COL))
    if var_type == "int":
        for i, iv in enumerate(field):
            fid.writelines(f"                                          {int(iv):d},  // {field_tag}[{i}]\n")
    elif var_type == "cusfloat":
        for i, iv in enumerate(field):
            fid.writelines(f"                                          {iv:0.16f},  // {field_tag}[{i}]\n")
    fid.writelines(f"                            " + "};\n")


def write_vector_header( fid, field: np.ndarray, field_tag: str, var_type: str, keyword="constexpr", is_inline=False )->None:
    inline_str = "inline" if is_inline else ""
    fid.writelines(str_start_col(f"alignas(FLOATING_PRECISION) {inline_str:s} static {keyword:s} {var_type}", f"{field_tag}[{field.shape[0]}] = " + "{\n", REF_COL))
    if var_type == "std::size_t":
        for i, iv in enumerate(field):
            fid.writelines(f"                                                                {int(iv):d},  // {field_tag}[{i}]\n")
    elif var_type == "cusfloat":
        for i, iv in enumerate(field):
            fid.writelines(f"                                                                {iv:0.16f},  // {field_tag}[{i}]\n")
    fid.writelines(f"                                                  " + "};\n")


def write_vector_line(fid, data: np.ndarray, field_name: str, var_type: str, int_name: str):
    if var_type == "cusfloat":
        field_data_str = ", ".join(f"{i:0.6E}" for i in data)
    elif var_type == "bool":
        field_data_str = ", ".join(f"{int(i):d}" for i in data)
    else:
        raise ValueError("Variable type not recognized!.")
    fid.writelines(str_start_col(f"{var_type}", f"{int_name}::{field_name}[{data.shape[0]}] = " + "{" + field_data_str + "};\n", REF_COL))
