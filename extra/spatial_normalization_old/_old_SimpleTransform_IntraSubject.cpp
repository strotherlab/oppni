#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
using namespace std;

// Updated 2013/07/29
//
// MakeTransform_IntraSubject: 
//
// generates transformation matrix to align EPI runs from the same subject
// e.g. rigid-body alignment
//
// Call as:   
//
// "g++ MakeTransform_IntraSubject.cpp"
//
// "./a.out {prefix for reference run} {prefix for run to be aligned}
// 
// ** Output produces a transformation matrix in the directory of aligned data,
//    with name < directory, "/spat_norm/Transmat_EPI_" , aligned-prefix,  "_to_" << reference-prefix << ".mat >

// ---------------------- main script -------------------------- //

int main(int argc, char* argv[])
{
    // defining key variables -----------------------------------
    int numvol, numsubj, idxX, idxY, returned;
    int ii, jj, kk, mm, nn;
	int idxIN1, idxIN2, idxOUT1, idxOUT2, idxSTRUC1, idxSTRUC2;
    //
    int iStart, iMiddl, iEnder, pp, scalNum, pptmp;
    int newSTART, newEND;
    //
    string refer_prefix, match_prefix, refer_sub1, refer_sub2, match_sub1, match_sub2;
    // ----------------------------------------------------------
    
     std::stringstream ss;
     std::string s;

    refer_prefix = argv[1]; // ARG#1: reference run prefix
    match_prefix = argv[2]; // ARG#2: matching run prefix
    
    string refer_string = string( refer_prefix );
    string match_string = string( match_prefix );
    
    // splitting
    idxX       = refer_string.find_last_of("/");
    refer_sub1 = string( refer_string, 0, idxX );
    refer_sub2 = string( refer_string,  idxX+1 );
    // splitting
    idxX       = match_string.find_last_of("/");
    match_sub1 = string( match_string, 0, idxX );
    match_sub2 = string( match_string,  idxX+1 );
    
    // .. create spatial-norming directory (in case it doesn't exist)
    ss.str(std::string());
    ss << "mkdir " << match_sub1 << "/spat_norm";
    s = ss.str();
    returned = system( s.c_str() ); 
    
    // .. compute temporal run-means, sent to spat-norm directory
    ss.str(std::string());
    ss << "fslmaths " << refer_string << "_baseproc.nii -Tmean " << refer_string << "_Tmean.nii";
    s = ss.str();
    returned = system( s.c_str() );                 
    //
    ss.str(std::string());
    ss << "fslmaths " << match_string << "_baseproc.nii -Tmean " << match_string << "_Tmean.nii";
    s = ss.str();
    returned = system( s.c_str() );                 
    // .. compute the spatial transform matrix + aligned run-mean, sent to aligned data directory
    ss.str(std::string());
    ss << "flirt -in " << refer_string << "_Tmean.nii -ref " << match_string << "_Tmean.nii -omat " << match_sub1 << "/spat_norm/Transmat_EPI_" << match_sub2 << "_to_" << refer_sub2 << ".mat -dof 6";
    s = ss.str();
    returned = system( s.c_str() );      
    
    // .. transform the mean, to check alignment
    ss.str(std::string());     
    ss<<"applyxfm4D " << match_string << "_Tmean.nii " << refer_string << "_Tmean.nii " << match_sub1 << "/spat_norm/" << match_sub2 << "_Tmean_intra_align.nii " << match_sub1 << "/spat_norm/Transmat_EPI_" << match_sub2 << "_to_" << refer_sub2 << ".mat -singlematrix";
    s = ss.str();
    returned = system( s.c_str() );   
    // .. rename and put reference mean in aligned data directory, also to check on alignment
    ss.str(std::string());     
    ss<<"mv " << refer_string << "_Tmean.nii.gz " << match_sub1 << "/spat_norm/" << match_sub2 << "_Tmean_intra_align_ref.nii.gz";
    s = ss.str();
    returned = system( s.c_str() );   
    // .. delete unaligned Tmean volume
    ss.str(std::string());     
    ss << "rm " << match_string << "_Tmean.nii.gz";
    s = ss.str();
    returned = system( s.c_str() );   

     
	return 0;
}
