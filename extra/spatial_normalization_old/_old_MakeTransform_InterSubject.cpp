#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
using namespace std;

// Updated 2013/08/13
//
// MakeTransform_InterSubject (VER. 3.0): 
// runs group-level spatial normalization across subjects
// based on anatomical template, after running Pipeline_Step1.
//
// Call as:   
//
// "g++ MakeTransform_InterSubject.cpp"
//
// "./a.out {subject inputs}.txt {reference T1 anatomical}.nii {optimization metric}"
// 
// where data specified in {subject inputs} are spatially matched to 
// {reference T1 anatomical}.nii template
// 
// {optimization metric}: abbreviation for metric used to perform pipeline 
// selection in STEP 3. e.g. for reproducibility, use "R"
// 
// ** Output results include spatial transformation matrix 
//         <OUTstr_sub1,"/spat_norm/Transmat_EPItoREF_", OUTstr_sub2n ".mat">
//
//    and spatially normed:
//    "base_proc" 4D epi data, with "..._sNorm.nii" appended to data name
//    optimized STEP3 pipeline SPMs
//    optimally preprocessed post-STEP3 4D data
//

// ------------------------------ functions ------------------------------ //

// FUNCTION: for "target" (dir/niivol.nii), get # volumes
int get_numvols( string file )
{
     int returned, ii, idxX, numvol;
     //
     std::stringstream ss;
     std::string s;
     std::string stringval;
     // 
     ss.str(std::string());
     ss << "nifti_tool -disp_nim -field nt -infiles " << file << " > temp_text.txt";
     s = ss.str();
     returned = system( s.c_str() );
     // ii) 
     ifstream temp ( "temp_text.txt" ); 
     while ( temp.good() )
     {
          // only update if line=good
          s = stringval;
          // query next line
          getline ( temp, stringval );
     }

     ii = 0;
     while( ii != s.npos ) { idxX = ii; ii = s.find(" ", idxX+1); }
     numvol = atoi( string(s, idxX).c_str() );
     
     return( numvol );
}

// ---------------------- main script -------------------------- //

int main(int argc, char* argv[])
{
    // defining key variables -----------------------------------
    int numsubj, idxX, idxY, returned;
	int idxOUT1, idxOUT2, idxSTRUC1, idxSTRUC2;
    //
    string stringoutline, sx;
    string intextname, inrefanatomic, metric_name;
    string OUTstr, STRUCstr;
    string OUTstr_sub1, OUTstr_sub2;
    string STRUCstr_sub1, STRUCstr_sub2;
    // ----------------------------------------------------------
    
    intextname    = argv[1];  // ARG#1: subject input textfile
    inrefanatomic = argv[2];  // ARG#2: name of T1 anatomical template
    metric_name   = argv[3];  // ARG#3: metric used to perform optimization
    
    // declare file stream:
    ifstream file_input ( intextname.c_str() ); 
    // declare the anatomical
    string refVolString = string( inrefanatomic );
    
    if( file_input.is_open()   ) // start reading in lines of input data
    {         
        //
        numsubj  = 0;

        while ( file_input.good() ) // parse textfile lines while something can be found
        {
            //// -------------------------------- ////
            ////   READING SUBJECTS INFORMATION   ////
            //// -------------------------------- ////

             getline ( file_input, stringoutline ); // read each line separately
             // parse input information:
             idxOUT1   = stringoutline.find("OUT="); 
             idxSTRUC1 = stringoutline.find("STRUCT=");   

             if ( (idxOUT1!=stringoutline.npos) & (idxSTRUC1!=stringoutline.npos) )
             {
                 numsubj++;
                 // 
                 // now increment to give starting index
                 idxOUT1  += 4;
                 idxSTRUC1+= 7; // for group norming
                 // index whitespace immediately following substring
                 idxOUT2  = stringoutline.find(" ",idxOUT1 +1);
                 idxSTRUC2= stringoutline.find(" ",idxSTRUC1 +1);  // for group norming
                 // 
                 OUTstr   = string( stringoutline, idxOUT1,  idxOUT2-idxOUT1 );
                 STRUCstr = string( stringoutline, idxSTRUC1,idxSTRUC2-idxSTRUC1 );  // for group norming
                 // splitting
                 idxX        = OUTstr.find_last_of("/");
                 OUTstr_sub1 = string( OUTstr, 0, idxX );
                 OUTstr_sub2 = string( OUTstr,  idxX+1 );
                 // splitting -- for group norming
                 idxX          = STRUCstr.find_last_of("/");
                 idxY          = STRUCstr.find_last_of(".");
                 STRUCstr_sub1 = string( STRUCstr,      0, idxX );
                 STRUCstr_sub2 = string( STRUCstr, idxX+1, idxY-idxX-1 );
                 // define string entities in the scope of this subject
                 std::stringstream ss;
                 std::string s;

                 // paths to data
                 ss.str(std::string());
                 ss <<OUTstr_sub1<<"/masks/"<<OUTstr_sub2<<"_mask.nii";
                 string maskpath_m1 = ss.str();

// ----------------------- GROUP NORMALIZATION BEGINS -----------------------  //

                 // get the "baseproc" volume as 4D fMRI reference
                 ss.str(std::string());
                 ss <<OUTstr<<"_baseproc";
                 string epiVoluString = ss.str();
                 string epiMaskString = maskpath_m1;                 
                 
                 // construct spatial normalization direcory
                 ss.str(std::string());
                 ss << "mkdir " <<OUTstr_sub1<<"/spat_norm";
                 s = ss.str();
                 returned = system( s.c_str() );

                 // begin by specifying some base parameters (e.g. voxel dimensions / identity transform)
                 // do only once, assuming that all subjects have same voxel dimensions, FOV etc.
                 //
                 if( numsubj == 1 )
                 {
                     cout << "running spatial-norm preparatory steps" << endl;

                     // output text file with acquisition parameters
                     ss.str(std::string());
                     ss << "fslhd " << epiVoluString << ".nii > niidims.txt"; 
                     s = ss.str();
                     returned = system( s.c_str() );

                     // initialize voxel dimension parameters
                     string dim1,dim2,dim3,dim4;
                     string pixdim1,pixdim2,pixdim3,pixdim4;
                     // open textfile
                     ifstream temp ( "niidims.txt" ); 
                     //
                     int kline    = 0;
                     int skipnext = 0;
                     // 
                     while ( temp.good() ) // run through the file
                     {

                         getline ( temp, sx );
                         int idx, idxSP;
                         kline++;

                         if( skipnext==0 ) // find lines with "dim<k>" --> #of voxels per dimension
                         {
                             idx = sx.find("dim1");   
                             if( idx!=s.npos ) { idxSP = sx.find_last_of(" "); dim1 = string( sx, idxSP );               }
                             idx = sx.find("dim2");   
                             if( idx!=sx.npos ) { idxSP = sx.find_last_of(" "); dim2 = string( sx, idxSP );              }
                             idx = sx.find("dim3");   
                             if( idx!=sx.npos ) { idxSP = sx.find_last_of(" "); dim3 = string( sx, idxSP );              }
                             idx = sx.find("dim4");   
                             if( idx!=sx.npos ) { idxSP = sx.find_last_of(" "); dim4 = string( sx, idxSP ); skipnext = 1;}
                         }
                         //
                         if( skipnext == 1 ) // find lines with "pixdim<k>" --> voxel widths (mm) / (sec on dim4)
                         {
                             idx = sx.find("pixdim1");   
                             if( idx!=sx.npos ) { idxSP = sx.find_last_of(" "); pixdim1 = string( sx, idxSP ); }
                             idx = sx.find("pixdim2");   
                             if( idx!=sx.npos ) { idxSP = sx.find_last_of(" "); pixdim2 = string( sx, idxSP ); }
                             idx = sx.find("pixdim3");   
                             if( idx!=sx.npos ) { idxSP = sx.find_last_of(" "); pixdim3 = string( sx, idxSP ); }
                             idx = sx.find("pixdim4");   
                             if( idx!=sx.npos ) { idxSP = sx.find_last_of(" "); pixdim4 = string( sx, idxSP ); }
                         }
                     }
                     temp.close();

                     // use information from "niidims.txt" to create a blank reference header
                     // this is used as common reference when transforming, so that all transformed volumes have same origin
                     ss.str(std::string());
                     ss << "fslcreatehd " << dim1<<" "<<dim2<<" "<<"35"<< " 1 " <<pixdim1<<" "<<pixdim2<<" "<<pixdim3<<" "<<pixdim4<<" 0 0 0 16 blankvol";
                     string s = ss.str();
                     returned = system( s.c_str() );
                     // make a unit matrix, used to do downsampling without transformation
                     returned = system( "rm eye.mat"  );
                     returned = system( "echo 1 0 0 0 >> eye.mat"  );
                     returned = system( "echo 0 1 0 0 >> eye.mat"  );
                     returned = system( "echo 0 0 1 0 >> eye.mat"  );
                     returned = system( "echo 0 0 0 1 >> eye.mat"  );
                 }

                 
                 // -------------------------[ brain-masking the T1 volume ]
                 ss.str(std::string());
                 ss << "3dSkullStrip -prefix " << OUTstr_sub1 << "/spat_norm/" << STRUCstr_sub2 << "_strip.nii -input " << STRUCstr;
                 s = ss.str();
                 returned = system( s.c_str() );       

                 cout << "Running normalization for :" << OUTstr << "..." << endl;

                 // -------------------------[ build transform --> T1 to Ref]
                 // transformation is 3D affine (translate/rotate/shear/reflect) with sinc interpolation
                 ss.str(std::string());
                 ss << "flirt -in " << OUTstr_sub1 << "/spat_norm/" << STRUCstr_sub2 << "_strip.nii -ref " << refVolString << " -out " << OUTstr_sub1 << "/spat_norm/" << STRUCstr_sub2 << "_T1toREF.nii -omat " << OUTstr_sub1 << "/spat_norm/Transmat_T1toREF_" << STRUCstr_sub2 << ".mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp sinc -sincwidth 7 -sincwindow hanning";
                 s = ss.str();
                 returned = system( s.c_str() );             
                 // -------------------------[ build transform --> EPI to T1]
                 // transformation is of EPI mask to masked T1 volume, using rigidbody+scaling (dof 7) transform
                 ss.str(std::string());
                 ss << "flirt -in " << epiMaskString << " -ref " << OUTstr_sub1 << "/spat_norm/" << STRUCstr_sub2 << "_strip.nii -omat " << OUTstr_sub1 << "/spat_norm/Transmat_EPItoT1_" << OUTstr_sub2 << ".mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 7  -interp sinc -sincwidth 7 -sincwindow hanning";
                 s = ss.str();
                 returned = system( s.c_str() );             
                 // -------------------------[ net matrix --> EPI to Ref ]
                 // concatenates (EPI to T1)x(T1 to Ref) matrices
                 ss.str(std::string());
                 ss << "convert_xfm -omat " << OUTstr_sub1 << "/spat_norm/Transmat_EPItoREF_" << OUTstr_sub2 << ".mat -concat " << OUTstr_sub1 << "/spat_norm/Transmat_T1toREF_" << STRUCstr_sub2 << ".mat " << OUTstr_sub1 << "/spat_norm/Transmat_EPItoT1_" << OUTstr_sub2 << ".mat";
                 s = ss.str();
                 returned = system( s.c_str() );             
                 // -------------------------[ build downsampled T1 for reference ]
                 ss.str(std::string());
                 ss << "flirt -in " << OUTstr_sub1 << "/spat_norm/" << STRUCstr_sub2 << "_T1toREF.nii -applyxfm -interp sinc -ref blankvol.nii.gz -init eye.mat -out " << OUTstr_sub1 << "/spat_norm/" << STRUCstr_sub2 << "_T1toREF_downsamp.nii";
                 s = ss.str();
                 returned = system( s.c_str() );             
                 
                 cout << "Normalizing baseproessing run for :" << OUTstr << "..." << endl;

                 // APPLY TRANSFORM-1: "base-processing" 4D volume
                 cout << epiVoluString<< endl;
                 ss.str(std::string());
                 ss<<"applyxfm4D "<<epiVoluString <<".nii " << OUTstr_sub1 << "/spat_norm/" << STRUCstr_sub2 << "_T1toREF_downsamp.nii " << epiVoluString << "_sNorm.nii " << OUTstr_sub1 << "/spat_norm/Transmat_EPItoREF_" << OUTstr_sub2 << ".mat -singlematrix";
                 s = ss.str();
                 returned = system( s.c_str() );   

                 // unzip normed volumes //
                 ss.str(std::string()); 
                 ss << "gunzip " << epiVoluString << "_sNorm.nii.gz";
                 s = ss.str();
                 returned = system( s.c_str() );                         
                 
                 // make reference mask
                 cout << epiVoluString<< endl;
                 ss.str(std::string());
                 ss<<"3dAutomask -prefix " << OUTstr_sub1 << "/spat_norm/" << OUTstr_sub2 << "_mask_sNorm.nii " << epiVoluString << "_sNorm.nii";
                 s = ss.str();
                 returned = system( s.c_str() );                    
                 
                 // APPLY TRANSFORM-2: optimized SPMs
                 // .string for 4D volume of optimized SPMs
                 ss.str(std::string());
                 ss << OUTstr_sub1 << "/matfiles/niftis_" << OUTstr_sub2 << "/Images_" << OUTstr_sub2 << "_opt_" << metric_name << "_Std_Fix_Ind";
                 epiVoluString = ss.str();
                 // .run transformation
                 ss.str(std::string());
                 ss<<"applyxfm4D "<< epiVoluString <<".nii " << OUTstr_sub1 << "/spat_norm/" << STRUCstr_sub2 << "_T1toREF_downsamp.nii " << epiVoluString << "_sNorm.nii " << OUTstr_sub1 << "/spat_norm/Transmat_EPItoREF_" << OUTstr_sub2 << ".mat -singlematrix";
                 s = ss.str();
                 returned = system( s.c_str() );   
                 
// CURRENTLY TURNED OFF, AS IT IS V. TIME CONSUMING
//
//                  // APPLY TRANSFORM-3: optimally preprocessed data
//                  // .sub-string for 4D volume of optimized pipelines
//                  ss.str(std::string());
//                  ss << OUTstr_sub1 << "/matfiles/niftis_" << OUTstr_sub2 << "/Preprocessed_opt_" << metric_name << "_" << OUTstr_sub2;
//                  epiVoluString = ss.str();
//                  // .run transformation STD
//                  ss.str(std::string());
//                  ss<<"applyxfm4D "<< epiVoluString <<"_STD.nii " << OUTstr_sub1 << "/spat_norm/" << STRUCstr_sub2 << "_T1toREF_downsamp.nii " << epiVoluString << "_STD_sNorm.nii " << OUTstr_sub1 << "/spat_norm/Transmat_EPItoREF_" << OUTstr_sub2 << ".mat -singlematrix";
//                  s = ss.str();
//                  returned = system( s.c_str() );   
//                  // .run transformation FIX
//                  ss.str(std::string());
//                  ss<<"applyxfm4D "<< epiVoluString <<"_FIX.nii " << OUTstr_sub1 << "/spat_norm/" << STRUCstr_sub2 << "_T1toREF_downsamp.nii " << epiVoluString << "_FIX_sNorm.nii " << OUTstr_sub1 << "/spat_norm/Transmat_EPItoREF_" << OUTstr_sub2 << ".mat -singlematrix";
//                  s = ss.str();
//                  returned = system( s.c_str() );   
//                  // .run transformation IND
//                  ss.str(std::string());
//                  ss<<"applyxfm4D "<< epiVoluString <<"_IND.nii " << OUTstr_sub1 << "/spat_norm/" << STRUCstr_sub2 << "_T1toREF_downsamp.nii " << epiVoluString << "_IND_sNorm.nii " << OUTstr_sub1 << "/spat_norm/Transmat_EPItoREF_" << OUTstr_sub2 << ".mat -singlematrix";
//                  s = ss.str();
//                  returned = system( s.c_str() );   
//                  
                 
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
