
function( print_src_list src_list )
    foreach( item IN LISTS ${src_list} )
        message( STATUS "Source file: ${item}")
    endforeach( )
endfunction( )