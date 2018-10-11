/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2016 University of Strathclyde and Authors ------
-------------------- Author: Francesco Torre -------------------------
-------------- e-mail: francesco.torre@strath.ac.uk ------------------
*/

#include "AstroUtils/smartastro_utils.h"

using namespace smartastro;

template < class T >
bool utilities::import_tetra_mesh(const std::string &file_name,
                                     std::vector< std::vector<T> > &mesh_vertices,
                                     std::vector< std::vector<unsigned int> > &mesh_faces)
{
    std::ifstream source(file_name.c_str());

	mesh_vertices.clear();
    mesh_faces.clear();

    std::string type;
    std::vector<T> coords(3);
    std::vector<unsigned int> verts(3);

    if(source.is_open())
    {
        while(source >> type)
        {
            if(type == "v")
            {
                source >> coords[0];
                source >> coords[1];
                source >> coords[2];
                mesh_vertices.push_back(coords);

                // std::cout << "\nv";
                // for(unsigned int index2 = 0; index2 < coords.size(); ++index2)
                // {
                //     std::cout << " " << coords[index2];
                // }
            }
            else if(type == "f")
            {
                source >> verts[0];
                source >> verts[1];
                source >> verts[2];
                mesh_faces.push_back(verts);

                // std::cout << "\nf";
                // for(unsigned int index2 = 0; index2 < verts.size(); ++index2)
                // {
                //     std::cout << " " << verts[index2];
                // }
            }
            else
            {
                std::cout << "\nError!" << std::endl;
                return false;
            }

            // std::cout << "\nPress Key To Continue" << std::endl;
            // std::cin.get();
        }
    }

    // std::cout << "\nmesh_vertices.size(): " << mesh_vertices.size() << std::endl;
    // std::cout << "v";
    // for(unsigned int index2 = 0; index2 < coords.size(); ++index2)
    // {
    //     std::cout << " " << mesh_vertices.front()[index2];
    // }
    // std::cout << std::endl;
    // std::cout << "v";
    // for(unsigned int index2 = 0; index2 < coords.size(); ++index2)
    // {
    //     std::cout << " " << mesh_vertices.back()[index2];
    // }
    // std::cout << std::endl;
    
    // std::cout << "\nmesh_faces.size(): " << mesh_faces.size() << std::endl;
    // std::cout << "f";
    // for(unsigned int index2 = 0; index2 < verts.size(); ++index2)
    // {
    //     std::cout << " " << mesh_faces.front()[index2];
    // }
    // std::cout << std::endl;
    // std::cout << "f";
    // for(unsigned int index2 = 0; index2 < verts.size(); ++index2)
    // {
    //     std::cout << " " << mesh_faces.back()[index2];
    // }
    // std::cout << std::endl;

    return true;
}


std::string utilities::vectord2string(const std::vector<double> vec, int precision)
{
    std::vector<double>::const_iterator it;
    std::stringstream s;
    for( it = vec.begin(); it != vec.end(); ++it )
    {
        if( it != vec.begin() )
            s << " ";

        s << std::setprecision(precision) << *it;
    }

    return s.str();
}

std::string utilities::exec_popen(const std::string &cmd)
{
    int RESULT_BUFFER = 128;
    char buffer[RESULT_BUFFER];
    std::string result = "";
    char * cmd_cstr = new char [cmd.length()+1];
    std::strcpy(cmd_cstr, cmd.c_str());
    FILE* pipe = popen(cmd_cstr, "r");
    if (!pipe)
    {
        throw std::runtime_error("popen() failed!");
    }
    try {
        while (!feof(pipe)) {
            if (fgets(buffer, RESULT_BUFFER, pipe) != NULL)
                result += buffer;
        }
    } catch (...) {
        pclose(pipe);
        throw;
    }
    pclose(pipe);
    return result;
}

template <typename I>
static I random_element(I begin, I end)
{
    const unsigned long n = std::distance(begin, end);
    const unsigned long divisor = (RAND_MAX + 1) / n;

    unsigned long k;
    do { k = std::rand() / divisor; } while (k >= n);

    std::cout << k << std::endl;
    std::advance(begin, k);
    return begin;
}




/**
 * Compare two strings neglecting case sensitive differences
 * @param str1: first string to compare
 * @param str2: second string to compare
 * @return true if equal, false otherwise
 */
bool smartastro::caseInsensitiveStringEqual( const std::string& str1, const std::string& str2 )
{
    // Define function to compare
    std::function<bool(char,char)> caseInsensitiveCharEqual =
            [](char a, char b) -> bool { return (std::toupper(a)==std::toupper(b)); };

    // Return bool
    return (
            (str1.size()==str2.size()) &&
            equal(str1.begin(),str1.end(),str2.begin(),caseInsensitiveCharEqual)
            );
}