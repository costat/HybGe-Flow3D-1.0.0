#include <boost/filesystem.hpp>

using namespace boost::filesystem;

bool
find_file( const path& ProblemPath, \
           const std::string& fileName, \
           path& MeshPath );
