# find eigen.so

string(TOUPPER eigen module_upper)
        execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c" 
                "import eigen; import os; print(eigen.__file__)"
                RESULT_VARIABLE _eigen_status 
                OUTPUT_VARIABLE _eigen_location
                ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
        if(NOT _eigen_status)
                set(PY_EIGEN_FILE_PATH ${_eigen_location} CACHE STRING 
                        "Location of Python module eigen")
        endif(NOT _eigen_status)
find_package_handle_standard_args(PY_eigen DEFAULT_MSG PY_EIGEN)