option(
  USE_CLANG_TIDY
  "Use clang-tidy if the program is found."
  OFF
)

if(USE_CLANG_TIDY)
  find_program(CLANG_TIDY clang-tidy)
  if(CLANG_TIDY)
    file(DOWNLOAD
      https://raw.githubusercontent.com/llvm-mirror/clang-tools-extra/master/clang-tidy/tool/run-clang-tidy.py
      ${PROJECT_BINARY_DIR}/run-clang-tidy.py
    )
    find_program(RUN_CLANG_TIDY run-clang-tidy.py)
    if(RUN_CLANG_TIDY)
      message("Using clang-tidy. Creating target... To run, use: make clang-tidy")
      add_custom_target(
        clang-tidy
        COMMAND python3 run-clang-tidy.py -quiet 2>&1
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        )
    else()
      message(
        FATAL_ERROR
        "run-clang-tidy.py not found. (Download and place in PATH). Aborting...")
    endif()
  else()
    message(FATAL_ERROR "clang-tidy not found. Aborting...")
  endif()
endif()
