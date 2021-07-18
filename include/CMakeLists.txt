# Install header files
install(DIRECTORY . DESTINATION include FILES_MATCHING PATTERN "*.hpp")

# Prepare and install Version.hpp
configure_file("pomerol/Version.hpp.in" "pomerol/Version.hpp")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/pomerol/Version.hpp"
        DESTINATION include/pomerol
)
