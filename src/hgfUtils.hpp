#include <vector>
#include <boost/filesystem.hpp>
#include "hgfAuxTools.hpp"

namespace bfs = boost::filesystem;

bool
find_file( const bfs::path& ProblemPath, \
           const std::string& fileName, \
           bfs::path& MeshPath );

void
geoFromFile( ProbParam& Par, const bfs::path& Geo );

void
problemParameters( ProbParam& Par, const bfs::path& ProblemPath);
