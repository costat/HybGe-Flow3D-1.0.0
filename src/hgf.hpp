#include <boost/filesystem.hpp>
#include "hgfAuxTools.hpp"

namespace bfs = boost::filesystem;

void
hgfDrive( const bfs::path& ProblemPath, const bfs::path& MeshPath, ProbParam& Par );
