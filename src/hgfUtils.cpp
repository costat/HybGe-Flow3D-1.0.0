#include <vector>
#include <iostream>
#include <boost/filesystem.hpp>

using namespace boost::filesystem;

/* function that checks if a file is present in a directory.
   used to determine if the mesh has already been built
   for the problem in this directory */
bool
find_file( const path& ProblemPath, \
           const std::string& fileName, \
           path& FilePath )
{
  if ( !exists( ProblemPath ) ) return false;
  directory_iterator end_itr;
  for ( directory_iterator itr( ProblemPath );
        itr != end_itr;
        ++itr ) {
    if ( is_directory(itr->status()) ) {
      if ( find_file(itr->path(), fileName, FilePath ) ) return true;
    }
    else if ( itr->path().leaf() == fileName ) {
      FilePath = itr->path();
      return true;
    }
  }
  return false;
}
