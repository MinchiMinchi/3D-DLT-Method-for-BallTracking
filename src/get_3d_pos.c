#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>


#define LIMIT_REF_POINTS 20
#define LIMIT_FLAME 6000
#define LIMIT_CAM 10


double REF_PixcelPos[LIMIT_CAM][LIMIT_REF_POINTS][2];
double REF_RealPos[LIMIT_REF_POINTS][3];


int NumOfSample;
int NumOfRef = 16;
int NumOfFrames = 0;

int NumOfCam = 2; // Num. of cameras


double PixcelPos[LIMIT_CAM][LIMIT_FLAME][2]; // ball pos. in 2d-image
double RealPos[LIMIT_FLAME][3]; // ball pos.


#define NUM_OF_DLT_PARAM 11 // DLT変数の数
double DLT_Param[LIMIT_CAM][NUM_OF_DLT_PARAM];


double DLT_x[NUM_OF_DLT_PARAM]; // DLT Parameters
double DLT_b[LIMIT_CAM][LIMIT_REF_POINTS*2];
double DLT_A[LIMIT_CAM][LIMIT_REF_POINTS*2][NUM_OF_DLT_PARAM];
double DLT_AtA[NUM_OF_DLT_PARAM][NUM_OF_DLT_PARAM];

double tmpVec[NUM_OF_DLT_PARAM];


double TransMat[LIMIT_CAM*2][3];
double TransVec[LIMIT_CAM*2];



void gauss_solve(double A[NUM_OF_DLT_PARAM][NUM_OF_DLT_PARAM],
		 double b[NUM_OF_DLT_PARAM],
		 int num)
{
  int i, j, k;
  double sigma;
    
  /* forward elimination */
  for(k=0 ; k<num-1 ; k++){

    if(A[k][k] < 1e-16){
      fprintf(stdout, "--- A[k][k] is Zero!! ----\n");
      fflush(stdout);
    }

    for(i=k+1 ; i<num ; i++){
      for(j=k+1 ; j<num ; j++){
	      A[i][j] = A[i][j]-A[k][j]*A[i][k]/A[k][k];
      }
      b[i] = b[i]-A[i][k]*b[k]/A[k][k];
    }
  }

  /* back substitution */
  b[num-1] = b[num-1] / A[num-1][num-1];
  for(i=num-2 ; i>=0 ; i--){
    sigma = 0.0;
    for(j=i+1 ; j<num ; j++){
      sigma += A[i][j]*b[j];
    }
    b[i] = 1/A[i][i]*(b[i]-sigma);
  }
}


void gauss_solve_3x3(double A[3][3],
		 double b[3])
{
  int i, j, k;
  double sigma;
    
  /* forward elimination */
  for(k=0 ; k<2 ; k++){

    if(A[k][k] < 1e-16){
      fprintf(stderr, "A[k][k] is Zero!!");
      fflush(stderr);
    }

    for(i=k+1 ; i<3 ; i++){
      for(j=k+1 ; j<3 ; j++){
      	A[i][j] = A[i][j]-A[k][j]*A[i][k]/A[k][k];
      }
      b[i] = b[i]-A[i][k]*b[k]/A[k][k];
    }
  }

  /* back substitution */
  b[2] = b[2] / A[2][2];
  for(i=1 ; i>=0 ; i--){
    sigma = 0.0;
    for(j=i+1 ; j<3 ; j++){
      sigma += A[i][j]*b[j];
    }
    b[i] = 1/A[i][i]*(b[i]-sigma);
  }
}


void calc_DLT_param()
{
  int i, j, k;
  int iCam;
  char str[256];
  FILE *fp;
  double k1, k2, b1, b2;
  char fname_cam[256] = "";
  char fname_3d[256] = "real-ref-points.csv";



// Read Kijun 3D
  fp = fopen(fname_3d, "r");
  if(fp == NULL){
    fprintf(stderr, "Can't Open File: %s\n", fname_3d);
    fflush(stderr);
    exit(1);
  }

  printf("----- 3D Kijun ------\n");
  fflush(stdout);

  /*
    fscanf(fp, "%s\n", str);
    printf("line1= %s\n", str);
    fflush(stdout);
  */
  for(i=0 ; i<NumOfRef ; i++){
    fscanf(fp, "%[^,],%lf,%lf,%lf", 
    str,
	  &REF_RealPos[i][0],
    &REF_RealPos[i][1],
    &REF_RealPos[i][2]);
  }

  for(i=0 ; i<NumOfRef ; i++){
    printf("P%d:  (%e, %e, %e) \n",
	   i, 
     REF_RealPos[i][0], 
     REF_RealPos[i][1],
     REF_RealPos[i][2]); 
  }

  fclose(fp);




  // Read Cam Ref. Points
  for(iCam=0 ; iCam<NumOfCam ; iCam++){
    sprintf(fname_cam, "cam%d-ref-points.csv", iCam+1);

    fp = fopen(fname_cam, "r");
    if(fp == NULL){
      fprintf(stderr, "Can't Open File: %s\n", fname_cam);
      fflush(stderr);
      exit(1);
    }

    printf("----- Cam %d ------\n", iCam);
    fflush(stdout);

    for(i=0 ; i<NumOfRef ; i++){
      fscanf(fp, "%[^,],%lf,%lf", 
	      str,
	      &REF_PixcelPos[iCam][i][0],
	      &REF_PixcelPos[iCam][i][1]);
    }


    for(i=0 ; i<NumOfRef ; i++){
      printf("P%d:  (%e, %e) \n",
	    i,
      REF_PixcelPos[iCam][i][0],
      REF_PixcelPos[iCam][i][1]);
    }

    fclose(fp);
  }



  // Calc.  Param.
  for(iCam=0 ; iCam<NumOfCam ; iCam++){
    for(i=0 ; i<NumOfRef ; i++){
      DLT_A[iCam][i*2][0] = REF_RealPos[i][0]; 
      DLT_A[iCam][i*2+1][0] = 0.0;

      DLT_A[iCam][i*2][1] = REF_RealPos[i][1]; ;
      DLT_A[iCam][i*2+1][1] = 0.0; 

      DLT_A[iCam][i*2][2] = REF_RealPos[i][2]; ;
      DLT_A[iCam][i*2+1][2] = 0.0;

      DLT_A[iCam][i*2][3] = 1.0;
      DLT_A[iCam][i*2+1][3] = 0.0;

      DLT_A[iCam][i*2][4] = 0.0;
      DLT_A[iCam][i*2+1][4] = REF_RealPos[i][0]; 

      DLT_A[iCam][i*2][5] = 0.0; 
      DLT_A[iCam][i*2+1][5] = REF_RealPos[i][1]; ;

      DLT_A[iCam][i*2][6] = 0.0;
      DLT_A[iCam][i*2+1][6] = REF_RealPos[i][2]; ;

      DLT_A[iCam][i*2][7] = 0.0;
      DLT_A[iCam][i*2+1][7] = 1.0; 

      DLT_A[iCam][i*2][8] = - REF_RealPos[i][0]*REF_PixcelPos[iCam][i][0]; 
      DLT_A[iCam][i*2+1][8] = - REF_RealPos[i][0]*REF_PixcelPos[iCam][i][1]; 

      DLT_A[iCam][i*2][9] = - REF_RealPos[i][1]*REF_PixcelPos[iCam][i][0]; 
      DLT_A[iCam][i*2+1][9] = - REF_RealPos[i][1]*REF_PixcelPos[iCam][i][1]; 

      DLT_A[iCam][i*2][10] = - REF_RealPos[i][2]*REF_PixcelPos[iCam][i][0]; 
      DLT_A[iCam][i*2+1][10] = - REF_RealPos[i][2]*REF_PixcelPos[iCam][i][1]; 
    
      DLT_b[iCam][i*2] = REF_PixcelPos[iCam][i][0];
      DLT_b[iCam][i*2+1] = REF_PixcelPos[iCam][i][1];
    }

    for(i=0 ; i<NUM_OF_DLT_PARAM ; i++){
      for(j=0 ; j<NUM_OF_DLT_PARAM ; j++){
        DLT_AtA[i][j] = 0.0;
        for(k=0 ; k<NumOfRef*2 ; k++){
	        DLT_AtA[i][j] += DLT_A[iCam][k][i]*DLT_A[iCam][k][j];
        }
      }
      DLT_Param[iCam][i] = 0.0;
      for(k=0 ; k<NumOfRef*2 ; k++){
        DLT_Param[iCam][i] += DLT_A[iCam][k][i]*DLT_b[iCam][k];
      }
    }

    // printf("---- pre calc. ----\n");
    // print_DLT_Param();

    for(i=0 ; i<NUM_OF_DLT_PARAM; i++){
      tmpVec[i] = DLT_Param[iCam][i];
    }

    // print_AtA(DLT_AtA, DLT_NUM_X);

    gauss_solve(DLT_AtA, tmpVec, NUM_OF_DLT_PARAM);

    for(i=0 ; i<NUM_OF_DLT_PARAM; i++){
      DLT_Param[iCam][i] = tmpVec[i];
    }

    // printf("---- post calc. ----\n");
    // print_DLT_Param();
  }

  // print_DLT_A(DLT_A[0], N_REF_POINTS*2, DLT_NUM_X);
  // print_DLT_b(DLT_b[0], N_REF_POINTS*2);
  // print_DLT_A(DLT_A[1], N_REF_POINTS*2, DLT_NUM_X);
  // print_DLT_b(DLT_b[1], N_REF_POINTS*2);

}



void PxPos_to_RealPos(
  double real_pos[3], 
  double px_pos[LIMIT_CAM][2])
{
  int i, j, k;
  double Mat[3][3];
  double vec[3];


  for(i=0 ; i<NumOfCam ; i++){
    TransMat[i*2][0] = DLT_Param[i][0] -  DLT_Param[i][8]*px_pos[i][0];
    TransMat[i*2][1] = DLT_Param[i][1] -  DLT_Param[i][9]*px_pos[i][0];
    TransMat[i*2][2] = DLT_Param[i][2] -  DLT_Param[i][10]*px_pos[i][0];

    TransMat[i*2+1][0] = DLT_Param[i][4] -  DLT_Param[i][8]*px_pos[i][1];
    TransMat[i*2+1][1] = DLT_Param[i][5] -  DLT_Param[i][9]*px_pos[i][1];
    TransMat[i*2+1][2] = DLT_Param[i][6] -  DLT_Param[i][10]*px_pos[i][1];

    TransVec[i*2] = px_pos[i][0] - DLT_Param[i][3];
    TransVec[i*2+1] = px_pos[i][1] - DLT_Param[i][7];
  }

  for(i=0 ; i<3 ; i++){
    for(j=0 ; j<3 ; j++){
      Mat[i][j] = 0.0;
      for(k=0 ; k<NumOfCam*2 ; k++){
	      Mat[i][j] += TransMat[k][i]*TransMat[k][j];
      }
    }
    vec[i] = 0.0;
    for(k=0 ; k<NumOfCam*2 ; k++){
      vec[i] += TransMat[k][i]*TransVec[k];
    }
  }
  
    gauss_solve_3x3(Mat, vec);

  for(i=0 ; i<3 ; i++){
    real_pos[i] = vec[i];
  }
}


int calc_3d_pos()
{
  int i;
  int iCam;
  char buf1[256], buf2[256], buf3[256];
  char str[256];
  FILE *fp_in[LIMIT_CAM];
  FILE *fp_out;
  FILE *fp_plot;
  char fname_in[256];
  char fname_out[256] = "out_3d_points.csv";
  char fname_plot[256] = "out_3d_plot.dat";
  double pixcel_pos[LIMIT_CAM][2];
  double real[3];


  // Open Cam Point File
  for(iCam=0 ; iCam<NumOfCam ; iCam++){
    sprintf(fname_in, "cam%d-ball-points.csv", iCam+1);

    fp_in[iCam] = fopen(fname_in, "r");
    if(fp_in[iCam] == NULL){
      fprintf(stderr, "Can't Open File: %s\n", fname_in);
      fflush(stderr);
      exit(1);
    }

    fscanf(fp_in[iCam], "%[^,],%[^,],%s", buf1, buf2, buf3);
  
    printf("------- open cam%d-ball-points --------\n", iCam+1);
    fflush(stdout);
  }
  


  printf("------- measure 3d points --------\n");

  fp_out = fopen(fname_out, "w");
  if(fp_out == NULL){
    fprintf(stderr, "Can't Open File \"%s\"\n", fname_out);
    fflush(stderr);
    exit(1);
  }

  fp_plot = fopen(fname_plot, "w");
  if(fp_out == NULL){
    fprintf(stderr, "Can't Open File \"%s\"\n", fname_plot);
    fflush(stderr);
    exit(1);
  }




  fprintf(fp_out, "frame_num,x,y,z\n"); 
  fprintf(fp_plot, "# x,  y,  z\n"); 


  NumOfFrames = 0;
  int isContinue = 1;
  int iframe;
  while( isContinue ){

    for(iCam=0 ; iCam<NumOfCam ; iCam++){
      if( fscanf(fp_in[iCam], "%d,%lf,%lf %[\n]", 
      &iframe, &pixcel_pos[iCam][0], &pixcel_pos[iCam][1]) == EOF){
        isContinue = 0;
      }
    }

    if(isContinue == 0){
      printf("NumOfFrame=%d, iframe=%d\n", NumOfFrames, iframe);
      fflush(stdout);
      break;
    }

    NumOfFrames++;
    printf("FrameNumber=%d,  ", iframe);
    for(iCam=0 ; iCam<NumOfCam ; iCam++){ 
      printf("cam%d=(%f, %f),  ",
      iCam+1 ,pixcel_pos[iCam][0], pixcel_pos[iCam][1]);
    }
    printf("\n");
    fflush(stdout);


    PxPos_to_RealPos(real,  pixcel_pos);

    fprintf(fp_out, "%d,%e,%e,%e\n", 
    iframe, real[0],real[1],real[2]);

    fprintf(fp_plot, "%e \t %e \t %e\n", 
    real[0],real[1],real[2]);


  }

  // close Cam Point File
  for(iCam=0 ; iCam<NumOfCam ; iCam++){
    fclose(fp_in[iCam]);
  }

  fclose(fp_out);
  fclose(fp_plot);
      
  
  
}


int main()
{
  calc_DLT_param();
  calc_3d_pos();
  return 0;
}
