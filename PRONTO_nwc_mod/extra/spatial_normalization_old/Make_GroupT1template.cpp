#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
using namespace std;

// Updated 2013/07/19
//
// Make_GroupT1template: 
// produces an iteratively-averaged group approximation of T1 brain map
//
// Call as:   
//
// "g++ Make_GroupT1template.cpp"
//
// "./a.out {subject inputs}.txt {reference T1 anatomical}.nii {prefix of new anatomical average} {number of iters}"
// 
// ** Output creates template: {prefix of new anatomical average}_iter{number of iters}.nii
//    that can be used in future spatial normalization

// ---------------------- main script -------------------------- //

int main(int argc, char* argv[])
{
    // defining key variables -----------------------------------
    int numsubj, idxX, idxY, returned;
    int iter, ii, N_iters;
	int idxOUT1, idxOUT2, idxSTRUC1, idxSTRUC2;
    //
    string stringoutline;
    string intextname, inrefanatomic, inprefix, avgrefprefix, iterstring;
    string OUTstr, STRUCstr;
    string OUTstr_sub1, OUTstr_sub2;
    string STRUCstr_sub1, STRUCstr_sub2;
    string OUTstr_array[100];
    string STRUCstr_array[100];
    // ----------------------------------------------------------
    std::stringstream setlist;
    
    intextname    = argv[1];  // ARG#1: subject input textfile
    inrefanatomic = argv[2];  // ARG#2: name of T1 anatomical template
    avgrefprefix  = argv[3];  // ARG#3: prefix of new T1 anatomical average
    iterstring    = argv[4];  // ARG#4: number of iterations
    
    N_iters = atoi( iterstring.c_str() );
    
    cout << "Generating iterative group mask, over " << N_iters << " iterations " << endl;
    
    // now, open these files 
    idxX     = intextname.find(".txt");
    inprefix = string( intextname, 0, idxX );
    // declare file stream:
    ifstream file_input ( intextname.c_str() ); 
    // declare the anatomical
    string refVolString = string( inrefanatomic );
    string avgRegPrefix = string( avgrefprefix );
    
    //// (1) First-run matching
    
    if( file_input.is_open()   ) // start reading in lines of input data
    {         
        //
        numsubj  = 0;
        // reformat list for averaging
        setlist.str(std::string());

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
                 // splitting the string
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
                 // concat + store for later iterations
                 OUTstr_array[numsubj-1]   = OUTstr;
                 STRUCstr_array[numsubj-1] = STRUCstr;
                 
// ----------------------- ZEROTH NORMALIZATION BEGINS -----------------------  //

                 // construct spatial normalization direcory
                 ss.str(std::string());
                 ss << "mkdir " <<OUTstr_sub1<<"/spat_norm";
                 s = ss.str();
                 returned = system( s.c_str() );

                 cout << "Running normalization for :" << OUTstr_sub2 << "..." << endl;
                 // -------------------------[ brain-masking the T1 volume ]
                 ss.str(std::string());
                 ss << "3dSkullStrip -prefix " << OUTstr_sub1 << "/spat_norm/" << STRUCstr_sub2 << "_strip.nii -input " << STRUCstr;
                 s = ss.str();
                 returned = system( s.c_str() );                 
                 // -------------------------[ build transform --> T1 to Ref]
                 // transformation is 3D affine (translate/rotate/shear/reflect) with sinc interpolation
                 ss.str(std::string());
                 ss << "flirt -in " << OUTstr_sub1 << "/spat_norm/" << STRUCstr_sub2 << "_strip.nii -ref " << refVolString << " -out " << OUTstr_sub1 << "/spat_norm/" << STRUCstr_sub2 << "_T1toREF_iter0.nii -omat " << OUTstr_sub1 << "/spat_norm/Transmat_T1toREF_" << STRUCstr_sub2 << "_iter0.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp sinc -sincwidth 7 -sincwindow hanning";
                 s = ss.str();
                 returned = system( s.c_str() );             
                     
                 // concat the list
                 if( numsubj == 1 ){  setlist << "fslmerge -t temp_cat_4D.nii " << OUTstr_sub1 << "/spat_norm/" << STRUCstr_sub2 << "_T1toREF_iter0.nii";  }
                 else              {  setlist << " " << OUTstr_sub1 << "/spat_norm/" << STRUCstr_sub2 << "_T1toREF_iter0.nii";               }
                 
// ------------------- ZEROTH NORMALIZATION ENDS -------------------  //
                 
             }
        }
    }
    
    for( iter=1; iter<=N_iters; iter++ )
    {
     
        // define string entities in the scope of this subject
        std::stringstream ss;
        std::string s;
        
        // compute the average --> this is the reference template for iteration #<iter>
        s = setlist.str();
        returned = system( s.c_str() );
        //
        ss.str(std::string());
        ss << "fslmaths temp_cat_4D.nii.gz -Tmean " << avgRegPrefix << "_iter"<< iter <<".nii";
        s = ss.str();
        returned = system( s.c_str() );      
        //
        ss.str(std::string());
        ss << "rm temp_cat_4D.nii.gz";
        s = ss.str();
        returned = system( s.c_str() );             
        
        if(iter>1)
        {
            // delete the previous iteration of the mean anatomical
            ss.str(std::string());
            ss << "rm " << avgRegPrefix << "_iter"<< (iter-1) <<".nii.gz";
            s = ss.str();
            returned = system( s.c_str() );                  
        }
        
        // reformat list for averaging
        setlist.str(std::string());

        for( ii=0; ii<numsubj; ii++ )
        {
            
// ------------------- ITERATED NORMALIZATION BEGINS -------------------  //
            
             // concat + store for next iterations
             OUTstr = OUTstr_array[ii];
             STRUCstr = STRUCstr_array[ii];

             // splitting the string
             idxX        = OUTstr.find_last_of("/");
             OUTstr_sub1 = string( OUTstr, 0, idxX );
             OUTstr_sub2 = string( OUTstr,  idxX+1 );
             // splitting -- for group norming
             idxX          = STRUCstr.find_last_of("/");
             idxY          = STRUCstr.find_last_of(".");
             STRUCstr_sub1 = string( STRUCstr,      0, idxX );
             STRUCstr_sub2 = string( STRUCstr, idxX+1, idxY-idxX-1 );

             // -------------------------[ build transform --> T1 to Ref]
             // transformation is 3D affine (translate/rotate/shear/reflect) with sinc interpolation
             ss.str(std::string());
             ss << "flirt -in " << OUTstr_sub1 << "/spat_norm/" << STRUCstr_sub2 << "_strip.nii -ref " << avgRegPrefix << "_iter"<< iter <<".nii -out " << OUTstr_sub1 << "/spat_norm/" << STRUCstr_sub2 << "_T1toREF_iter"<< iter <<".nii -omat " << OUTstr_sub1 << "/spat_norm/Transmat_T1toREF_" << STRUCstr_sub2 << "_iter"<< iter <<".mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp sinc -sincwidth 7 -sincwindow hanning";
             s = ss.str();
             returned = system( s.c_str() );     
             
             // delete the previous iteration of the transformed anatomicals
             // This includes "zeroth" step outputs!
             ss.str(std::string());
             ss << "rm " << OUTstr_sub1 << "/spat_norm/" << STRUCstr_sub2 << "_T1toREF_iter"<< (iter-1) <<".nii.gz";
             s = ss.str();
             returned = system( s.c_str() );                     
             //.... and the transformation matrix
             ss.str(std::string());
             ss << "rm " << OUTstr_sub1 << "/spat_norm/Transmat_T1toREF_" << STRUCstr_sub2 << "_iter"<< (iter-1) <<".mat";
             s = ss.str();
             returned = system( s.c_str() );    

             // concat the list
             if( ii == 0 ){  setlist << "fslmerge -t temp_cat_4D.nii " << OUTstr_sub1 << "/spat_norm/" << STRUCstr_sub2 << "_T1toREF_iter"<< iter <<".nii";  }
             else         {  setlist << " " << OUTstr_sub1 << "/spat_norm/" << STRUCstr_sub2 << "_T1toREF_iter"<< iter <<".nii";               }

// ------------------- ITERATED NORMALIZATION ENDS -------------------  //
        }
    }

 
	return 0;
}
