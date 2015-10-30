// Open the mesh named in command line arguments,
// update elems in oldid to newid

#include <iostream>

#include "libmesh.h"

#include "mesh.h"
#include "elem.h"
#include "getpot.h"
#include "boundary_info.h"

void usage_error(const char *progname)
{
  std::cout << "Usage: " << progname
            << " --input inputmesh --output outputmesh --oldid --newid" 
            << std::endl;

  exit(1);
}

int main(int argc, char** argv)
{
  LibMeshInit init(argc, argv);

  GetPot cl(argc, argv);

  Mesh mesh;

  if(!cl.search("--input"))
  {
    std::cerr << "No --input argument found!" << std::endl;
    usage_error(argv[0]);
  } 
  const char* meshname = cl.next("");

  mesh.read(meshname);
  std::cout << "Loaded mesh " << meshname << std::endl;

  if(!cl.search("--oldid"))
  {
    std::cerr << "No --oldid argument found!" << std::endl;
    usage_error(argv[0]);
  }
  unsigned int oldid = cl.next(0);

  if(!cl.search("--newid"))
  {
    std::cerr << "No --newid argument found!" << std::endl;
    usage_error(argv[0]);
  }
  unsigned int newid = cl.next(0);

  MeshBase::element_iterator           el = mesh.elements_begin();
  const MeshBase::element_iterator end_el = mesh.elements_end();
  
  for (; el != end_el; ++el)
  {
    Elem *elem = *el;

    // Update the elements in old subdomain to the new subdomain
    if (elem->subdomain_id() == oldid)
      elem->subdomain_id() = newid;

    // Update sides that match this id
    unsigned int n_sides = elem->n_sides();
    for (unsigned int s=0; s != n_sides; ++s)
    {
      std::vector<short int> boundary_ids = mesh.boundary_info->boundary_ids(elem, s);
      if (std::find(boundary_ids.begin(), boundary_ids.end(), oldid) != boundary_ids.end())
      {
        mesh.boundary_info->remove_side(elem, s, oldid);
        mesh.boundary_info->add_side(elem, s, newid);
      }
    }
  }
  // Finally update all the nodesets
  mesh.boundary_info->clear_boundary_node_ids();
  mesh.boundary_info->build_node_list_from_side_list();
//  mesh.prepare_for_use();

  std::string outputname;
  if(cl.search("--output"))
  {
    outputname = cl.next("");
  } 
  else
  { 
    outputname = "new.";
    outputname += meshname;
  } 

  mesh.write(outputname.c_str());
  std::cout << "Wrote mesh " << outputname << std::endl;

  return 0;
}
