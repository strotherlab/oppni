#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
using namespace std;

// Updated 2013/08/13
//
// MakeTransform_IntraSubject (VER. 2.0): 
// runs rigid-body spatial normalization, matching data runs in
// one textfile to paired runs in another
//
// Call as:   
//
// "g++ MakeTransform_IntraSubject.cpp"
//
// "./a.out {reference inputs}.txt {alignment inputs}.txt {optimization metric}"
// 
// where data specified in {alignment inputs} spatially matched to
// data in {reference inputs}. Each "alignment" line is matched to
// its corresponding line in "reference".
// e.g. data specified in line-1 (alignment) matched to line-1 (reference), etc.
// 
// {optimization metric}: abbreviation for metric used to perform pipeline 
// selection in STEP 3. e.g. for reproducibility, use "R"
// 
// ** Output results include spatial transformation matrix 
// <OUTstr_aln_sub1, "/spat_norm/Transmat_EPI_", OUTstr_aln_sub2, "_to_", OUTstr_ref_sub2, ".mat>
//
//    and spatially normed:
//      optimized STEP3 pipeline SPMs
//      optimally preprocessed post-STEP3 4D data
//

// ---------------------- main script -------------------------- //

int main(int argc, char* argv[])
{
    // defining key variables -----------------------------------
    int numsubj, idxX, returned;
	int idxOUT_aln1, idxOUT_aln2, idxOUT_ref1, idxOUT_ref2;
    //
    string stringout_ref, stringout_aln, sx;
    string reftextname, aligntextname, metric_name;
    string OUTstr_ref, OUTstr_aln;
    string OUTstr_ref_sub1, OUTstr_ref_sub2, OUTstr_aln_sub1, OUTstr_aln_sub2;
    // ----------------------------------------------------------
    
    reftextname   = argv[1];  // ARG#1: reference input textfile
    aligntextname = argv[2];  // ARG#2: align input textfile
    metric_name   = argv[3];  // ARG#3: metric used to perform optimization
    
    // declare file stream:
    ifstream file_ref (   reftextname.c_str() ); 
    ifstream file_aln ( aligntextname.c_str() ); 
    
    if( file_ref.is_open() & file_aln.is_open()   ) // start reading in lines of input data
    {         
        //
        numsubj  = 0;

        while ( file_ref.good() & file_aln.good() ) // parse textfile lines while something can be found
        {
            //// -------------------------------- ////
            ////   READING SUBJECTS INFORMATION   ////
            //// -------------------------------- ////

             getline ( file_ref, stringout_ref ); // read each line separately
             getline ( file_aln, stringout_aln ); // read each line separately
             // parse input information:
             idxOUT_ref1   = stringout_ref.find("OUT="); 
             idxOUT_aln1   = stringout_aln.find("OUT="); 
             
             if ( (idxOUT_ref1!=stringout_ref.npos) & (idxOUT_aln1!=stringout_aln.npos) )
             {
                 numsubj++;
                 // 
                 // now increment to give starting index
                 idxOUT_ref1  += 4;
                 idxOUT_aln1  += 4;
                 // index whitespace immediately following substring
                 idxOUT_ref2  = stringout_ref.find(" ",idxOUT_ref1 +1);
                 idxOUT_aln2  = stringout_ref.find(" ",idxOUT_aln1 +1);
                 // 
                 OUTstr_ref   = string( stringout_ref, idxOUT_ref1,  idxOUT_ref2-idxOUT_ref1 );
                 OUTstr_aln   = string( stringout_aln, idxOUT_aln1,  idxOUT_aln2-idxOUT_aln1 );
                 // splitting
                 idxX        = OUTstr_ref.find_last_of("/");
                 OUTstr_ref_sub1 = string( OUTstr_ref, 0, idxX );
                 OUTstr_ref_sub2 = string( OUTstr_ref,  idxX+1 );
                 //
                 idxX        = OUTstr_aln.find_last_of("/");
                 OUTstr_aln_sub1 = string( OUTstr_aln, 0, idxX );
                 OUTstr_aln_sub2 = string( OUTstr_aln,  idxX+1 );                 
                 // define string entities in the scope of this subject
                 std::stringstream ss;
                 std::string s;

// ----------------------- GROUP NORMALIZATION BEGINS -----------------------  //

                 // get the "baseproc" volume as 4D fMRI reference
                 ss.str(std::string());
                 ss <<OUTstr_aln<<"_baseproc";
                 string epiVoluString = ss.str();
                 
                 // construct spatial normalization direcory
                 ss.str(std::string());
                 ss << "mkdir " <<OUTstr_aln_sub1<<"/spat_norm";
                 s = ss.str();
                 returned = system( s.c_str() );

                //------------------------------------
                // building the rigid-body transforms
                 
                // .. compute temporal run-means, sent to spat-norm directory
                ss.str(std::string());
                ss << "fslmaths " << OUTstr_ref << "_baseproc.nii -Tmean " << OUTstr_ref << "_Tmean.nii";
                s = ss.str();
                returned = system( s.c_str() );                 
                //
                ss.str(std::string());
                ss << "fslmaths " << OUTstr_aln << "_baseproc.nii -Tmean " << OUTstr_aln << "_Tmean.nii";
                s = ss.str();
                returned = system( s.c_str() );                 
                // .. compute the spatial transform matrix + aligned run-mean, sent to aligned data directory
                ss.str(std::string());
                ss << "flirt -in " << OUTstr_ref << "_Tmean.nii -ref " << OUTstr_aln << "_Tmean.nii -omat " << OUTstr_aln_sub1 << "/spat_norm/Transmat_EPI_" << OUTstr_aln_sub2 << "_to_" << OUTstr_ref_sub2 << ".mat -dof 6";
                s = ss.str();
                returned = system( s.c_str() );      

                // .. transform the mean, to check alignment
                ss.str(std::string());     
                ss<<"applyxfm4D " << OUTstr_aln << "_Tmean.nii " << OUTstr_ref << "_Tmean.nii " << OUTstr_aln_sub1 << "/spat_norm/" << OUTstr_aln_sub2 << "_Tmean_intra_align.nii " << OUTstr_aln_sub1 << "/spat_norm/Transmat_EPI_" << OUTstr_aln_sub2 << "_to_" << OUTstr_ref_sub2 << ".mat -singlematrix";
                s = ss.str();
                returned = system( s.c_str() );   
                // .. rename and put reference mean in aligned data directory, also to check on alignment
                ss.str(std::string());     
                ss<<"mv " << OUTstr_ref << "_Tmean.nii.gz " << OUTstr_aln_sub1 << "/spat_norm/" << OUTstr_aln_sub2 << "_Tmean_intra_align_ref.nii.gz";
                s = ss.str();
                returned = system( s.c_str() );   
                // .. delete unaligned Tmean volume
                ss.str(std::string());     
                ss << "rm " << OUTstr_aln << "_Tmean.nii.gz";
                s = ss.str();
                returned = system( s.c_str() );   


                 // APPLY TRANSFORM-2: optimized SPMs
                 // .string for 4D volume of optimized SPMs
                 ss.str(std::string());
                 ss << OUTstr_aln_sub1 << "/matfiles/niftis_" << OUTstr_aln_sub2 << "/Images_" << OUTstr_aln_sub2 << "_opt_" << metric_name << "_Std_Fix_Ind";
                 epiVoluString = ss.str();
                 // .run transformation
                 ss.str(std::string());
                 ss<<"applyxfm4D " <<epiVoluString <<".nii " << OUTstr_aln_sub1 << "/spat_norm/" << OUTstr_aln_sub2 << "_Tmean_intra_align_ref.nii.gz " << epiVoluString << "_sAlign.nii " << OUTstr_aln_sub1 << "/spat_norm/Transmat_EPI_" << OUTstr_aln_sub2 << "_to_" << OUTstr_ref_sub2 << ".mat -singlematrix";
                 s = ss.str();
                 returned = system( s.c_str() );   
                 
                 // APPLY TRANSFORM-3: optimally preprocessed data
                 // .sub-string for 4D volume of optimized pipelines
                 ss.str(std::string());
                 ss << OUTstr_aln_sub1 << "/matfiles/niftis_" << OUTstr_aln_sub2 << "/Preprocessed_opt_" << metric_name << "_" << OUTstr_aln_sub2;
                 epiVoluString = ss.str();
                 // .run transformation STD
                 ss.str(std::string());
                 ss<<"applyxfm4D " <<epiVoluString <<"_STD.nii " << OUTstr_aln_sub1 << "/spat_norm/" << OUTstr_aln_sub2 << "_Tmean_intra_align_ref.nii.gz " << epiVoluString << "_STD_sAlign.nii " << OUTstr_aln_sub1 << "/spat_norm/Transmat_EPI_" << OUTstr_aln_sub2 << "_to_" << OUTstr_ref_sub2 << ".mat -singlematrix";
                 s = ss.str();
                 returned = system( s.c_str() );   
                 // .run transformation FIX
                 ss.str(std::string());
                 ss<<"applyxfm4D " <<epiVoluString <<"_FIX.nii " << OUTstr_aln_sub1 << "/spat_norm/" << OUTstr_aln_sub2 << "_Tmean_intra_align_ref.nii.gz " << epiVoluString << "_FIX_sAlign.nii " << OUTstr_aln_sub1 << "/spat_norm/Transmat_EPI_" << OUTstr_aln_sub2 << "_to_" << OUTstr_ref_sub2 << ".mat -singlematrix";
                 s = ss.str();
                 returned = system( s.c_str() );   
                 // .run transformation IND
                 ss.str(std::string());
                 ss<<"applyxfm4D " <<epiVoluString <<"_IND.nii " << OUTstr_aln_sub1 << "/spat_norm/" << OUTstr_aln_sub2 << "_Tmean_intra_align_ref.nii.gz " << epiVoluString << "_IND_sAlign.nii " << OUTstr_aln_sub1 << "/spat_norm/Transmat_EPI_" << OUTstr_aln_sub2 << "_to_" << OUTstr_ref_sub2 << ".mat -singlematrix";
                 s = ss.str();
                 returned = system( s.c_str() );   
                 
                 
// ------------------- GROUP NORMALIZATION ENDS -------------------  //

             }
        }            
    }
    else
    {
        cout << "ERROR: Unable to read/open input subjectlist - terminating early!\n" << endl;            
    }

 
	return 0;
}
