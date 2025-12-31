
function( generate_coeffs_header FREQ_DB_PREC )
    set( include_lines "" )
    list( APPEND include_lines "#include \"./fin_depth_coeffs/${FREQ_DB_PREC}/L0.hpp\""       )
    list( APPEND include_lines "#include \"./fin_depth_coeffs/${FREQ_DB_PREC}/L0_dA.hpp\""    )
    list( APPEND include_lines "#include \"./fin_depth_coeffs/${FREQ_DB_PREC}/L0_dB.hpp\""    )
    list( APPEND include_lines "#include \"./fin_depth_coeffs/${FREQ_DB_PREC}/L1.hpp\""       )
    list( APPEND include_lines "#include \"./fin_depth_coeffs/${FREQ_DB_PREC}/L1_dA.hpp\""    )
    list( APPEND include_lines "#include \"./fin_depth_coeffs/${FREQ_DB_PREC}/L1_dB.hpp\""    )
    list( APPEND include_lines "#include \"./fin_depth_coeffs/${FREQ_DB_PREC}/L2.hpp\""       )
    list( APPEND include_lines "#include \"./fin_depth_coeffs/${FREQ_DB_PREC}/L3.hpp\""       )
    list( APPEND include_lines "#include \"./fin_depth_coeffs/${FREQ_DB_PREC}/L3_dA.hpp\""    )
    list( APPEND include_lines "#include \"./fin_depth_coeffs/${FREQ_DB_PREC}/L3_dB.hpp\""    )
    list( APPEND include_lines "#include \"./fin_depth_coeffs/${FREQ_DB_PREC}/Linf.hpp\""     )
    list( APPEND include_lines "#include \"./fin_depth_coeffs/${FREQ_DB_PREC}/Linf_dA.hpp\""  )
    list( APPEND include_lines "#include \"./fin_depth_coeffs/${FREQ_DB_PREC}/Linf_dB.hpp\""  )
    list( APPEND include_lines "#include \"./fin_depth_coeffs/${FREQ_DB_PREC}/M1.hpp\""       )
    list( APPEND include_lines "#include \"./fin_depth_coeffs/${FREQ_DB_PREC}/M1_dA.hpp\""    )
    list( APPEND include_lines "#include \"./fin_depth_coeffs/${FREQ_DB_PREC}/M1_dB.hpp\""    )
    list( APPEND include_lines "#include \"./fin_depth_coeffs/${FREQ_DB_PREC}/M2.hpp\""       )
    list( APPEND include_lines "#include \"./fin_depth_coeffs/${FREQ_DB_PREC}/M3.hpp\""       )
    list( APPEND include_lines "#include \"./fin_depth_coeffs/${FREQ_DB_PREC}/M3_dA.hpp\""    )
    list( APPEND include_lines "#include \"./fin_depth_coeffs/${FREQ_DB_PREC}/M3_dB.hpp\""    )
    list( APPEND include_lines "#include \"./inf_depth_coeffs/${FREQ_DB_PREC}/R00_dX.hpp\""   )
    list( APPEND include_lines "#include \"./inf_depth_coeffs/${FREQ_DB_PREC}/R11.hpp\""      )
    list( APPEND include_lines "#include \"./inf_depth_coeffs/${FREQ_DB_PREC}/R11_dX.hpp\""   )

    string(JOIN "\n" include_text ${include_lines})

    file(
            WRITE 
            "${CMAKE_CURRENT_SOURCE_DIR}/src/green/fd_coeffs_files.hpp"
            "#pragma once\n\n"
            "// Auto-generated. Do not edit.\n\n"
            "${include_text}\n"
        )
    
endfunction()