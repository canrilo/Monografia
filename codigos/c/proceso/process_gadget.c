//#define _FILE_OFFSET_BITS 64
//#include <features.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw.h>
struct io_header_1
{
  int npart[6];
  double mass[6];
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  int npartTotal[6];
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8];	/* fills to 256 Bytes */
} header1;
int NumPart, Ngas, max=0;
struct particle_data
{
  float Pos[3];
  float Vel[3];
  float Mass;
  int Type;

  float Rho, U, Temp, Ne;
} *P;
int *Id,*G;
double Time, Redshift;
struct halo{
	int num;
	float posCM[3];
	float velCM[3];
}*H;
float masa=1;
char nombrefof[20]="fof_data_snp.dat";
char nombregrupos[10]="fof.grp";
float *densidades;

/*=================================================================================*/
void datos_halo(int ind);
int load_snapshot(char *fname, int files);
int allocate_memory(void);
int reordering(void);
void exportar_fof(void);
void centros(char *archivo);
void cic(int div);
void add_dens(int x, int y, int z, float peso, int div);
int sign(int v);
void escribir_densidades(char nom_dens[], int div);
/*=================================================================================*/

int main(void){
	char snapshot[100],arch_densidades[100];
	char comando[100], nombre_gr_final[70];
	int des,divisiones;
	char ll[10];
	char longi[10];
	float link;
	printf("Welcome to the Gadget-2 snapshot analysis program\n");
	printf("Please give the name of the snapshot to be used\n");
	scanf("%s",snapshot);
	load_snapshot(snapshot,1);
	reordering();
	printf("The number of particles is %d\n",NumPart);
	printf("The mass is %f\n",header1.mass[1]);
	printf("The size of the box is %.1f\n",header1.BoxSize);
	printf("Snapshot loaded. What would you like to do\n");
	printf("1. Power spectrum\n2. Halo Finder\n3. Export File\n------------------\n");
	scanf("%d",&des);
	switch (des)
	{
		case 1:
			printf("Enter the number of divisions for CIC\n");
			scanf("%d",&divisiones);
			cic(divisiones);
			printf("Enter a name for the file of densities:\n");
			scanf("%s",arch_densidades);
			escribir_densidades(arch_densidades,divisiones);
			break;
		case 2: 
			exportar_fof();
			link=header1.BoxSize*0.17/pow(NumPart,1/3.);
			printf("Enter a value for the linking lenght (Recomended %.2f):\n",link);
			scanf("%s",ll);
			snprintf(longi,sizeof(longi),"%.1f",header1.BoxSize);
			strcpy(comando,"./fof -e ");
			strcat(comando,ll);
			strcat(comando," -m 10 -px ");
			strcat(comando,longi);
			strcat(comando," -py ");
			strcat(comando,longi);
			strcat(comando," -pz ");
			strcat(comando,longi);
			strcat(comando," < ");
			strcat(comando,nombrefof);
			//printf("%s",comando);
			system(comando);
			printf("Enter a name for the group halo data file:\n");
			scanf("%s",nombre_gr_final);
			centros(nombre_gr_final);
			break;
		case 3:
			printf("Enter a filename\n");
			break;
		default:
			printf("You must give a valid number\n");
			
	}
	
	
}

/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.
 */
int load_snapshot(char *fname, int files)
{
  FILE *fd;
  char buf[200];
  int i, j, k, dummy, ntot_withmasses;
  int t, n, off, pc, pc_new, pc_sph;

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

  for(i = 0, pc = 1; i < files; i++, pc = pc_new)
    {
      if(files > 1)
	sprintf(buf, "%s.%d", fname, i);
      else
	sprintf(buf, "%s", fname);

      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open file `%s`\n", buf);
	  exit(0);
	}

      printf("reading `%s' ...\n", buf);
      fflush(stdout);
	
	  fread(&dummy, sizeof(dummy), 1, fd);
      fread(&header1, sizeof(header1), 1, fd);
      fread(&dummy, sizeof(dummy), 1, fd);
		
      if(files == 1)
		{
		  for(k = 0, NumPart = 0, ntot_withmasses = 0; k < 6; k++)
			NumPart += header1.npart[k];
		  Ngas = header1.npart[0];
		}
      else
		{
		  for(k = 0, NumPart = 0, ntot_withmasses = 0; k < 6; k++)
			NumPart += header1.npartTotal[k];
		  Ngas = header1.npartTotal[0];
		}

      for(k = 0, ntot_withmasses = 0; k < 6; k++)
	{
	  if(header1.mass[k] == 0)
	    ntot_withmasses += header1.npart[k];
	}

      if(i == 0)
	allocate_memory();

      SKIP;
      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header1.npart[k]; n++)
	    {
			//printf("%d\n",n);
	      fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
	      pc_new++;
	    }
	}
		
      SKIP;
		
      SKIP;
      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header1.npart[k]; n++)
	    {
	      fread(&P[pc_new].Vel[0], sizeof(float), 3, fd);
	      pc_new++;
	    }
	}
      SKIP;


      SKIP;
      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header1.npart[k]; n++)
	    {
	      fread(&Id[pc_new], sizeof(int), 1, fd);
	      pc_new++;
	    }
	}
      SKIP;


      if(ntot_withmasses > 0)
	SKIP;
      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header1.npart[k]; n++)
	    {
	      P[pc_new].Type = k;

	      if(header1.mass[k] == 0)
		fread(&P[pc_new].Mass, sizeof(float), 1, fd);
	      else
		P[pc_new].Mass = header1.mass[k];
	      pc_new++;
	    }
	}
      if(ntot_withmasses > 0)
	SKIP;


      if(header1.npart[0] > 0)
	{
	  SKIP;
	  for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
	    {
	      fread(&P[pc_sph].U, sizeof(float), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  SKIP;
	  for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
	    {
	      fread(&P[pc_sph].Rho, sizeof(float), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  if(header1.flag_cooling)
	    {
	      SKIP;
	      for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
		{
		  fread(&P[pc_sph].Ne, sizeof(float), 1, fd);
		  pc_sph++;
		}
	      SKIP;
	    }
	  else
	    for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
	      {
		P[pc_sph].Ne = 1.0;
		pc_sph++;
	      }
	}

      fclose(fd);
    }


  Time = header1.time;
  Redshift = header1.time;
}

/* this routine allocates the memory for the 
 * particle data.
 */
int allocate_memory(void)
{
  //printf("allocating memory...\n");

  if(!(P = malloc(NumPart * sizeof(struct particle_data))))
    {
      fprintf(stderr, "failed to allocate memory.\n");
      exit(0);
    }

  P--;				/* start with offset 1 */


  if(!(Id = malloc(NumPart * sizeof(int))))
    {
      fprintf(stderr, "failed to allocate memory.\n");
      exit(0);
    }

  Id--;				/* start with offset 1 */

  printf("allocating memory...done\n");
}

/* This routine brings the particles back into
 * the order of their ID's.
 * NOTE: The routine only works if the ID's cover
 * the range from 1 to NumPart !
 * In other cases, one has to use more general
 * sorting routines.
 */
int reordering(void)
{
  int i, j;
  int idsource, idsave, dest;
  struct particle_data psave, psource;


  printf("reordering....\n");

  for(i = 1; i <= NumPart; i++)
    {
      if(Id[i] != i)
	{
	  psource = P[i];
	  idsource = Id[i];
	  dest = Id[i];

	  do
	    {
	      psave = P[dest];
	      idsave = Id[dest];

	      P[dest] = psource;
	      Id[dest] = idsource;

	      if(dest == i)
		break;

	      psource = psave;
	      idsource = idsave;

	      dest = idsource;
	    }
	  while(1);
	}
    }

  printf("done.\n");

  Id++;
  free(Id);

  printf("space for particle ID freed\n");
}

/*Función que exporta el archivo con el formato para FOFhax*/
void exportar_fof(void){

	FILE *exp;
	int ind;
	
	exp=fopen(nombrefof,"w");
	if(!exp){
		printf("Hubo problemas abriendo el archivo %s\n",nombrefof);
	}
	printf("Escribiendo el archivo para FoF\n");
	fprintf(exp,"%d\n",NumPart);	 //Num total partículas
	fprintf(exp,"%d\n",NumPart);	//Num DM
	fprintf(exp,"0\n");				//Num gas
	fprintf(exp,"0\n");				//Num stars
	fprintf(exp,"0.01\n");			//time REVISARR!!!!!!!
	fprintf(exp,"0\n");				//nactive
	
	
	for(ind =1;ind<=NumPart;ind++){
		fprintf(exp,"%f %f %f\n",P[ind].Pos[0],P[ind].Pos[1],P[ind].Pos[2]);
	}
	printf("Archivo finalizado\n");
	fclose(exp);
}

/*Función que calcula centros de masa y velocidades con el índice de grupo*/
/*Se debe organizar para evitar hacer el recorrido tantas veces*/
void datos_halo(int ind){
	
	float x=0,y=0,z=0,vx=0,vy=0,vz=0;
	int j, cont=0;
	
	for(j=1;j<=NumPart;j++){
		if(G[j]==ind){
			x+=P[j].Pos[0];
			y+=P[j].Pos[1];
			z+=P[j].Pos[2];
			vx+=P[j].Vel[0];
			vy+=P[j].Vel[1];
			vz+=P[j].Vel[2];
			cont++;
		}
	}
	H[ind].posCM[0]=(float)x/cont;
	H[ind].posCM[1]=(float)y/cont;
	H[ind].posCM[2]=(float)z/cont;
	H[ind].velCM[0]=(float)vx/cont;
	H[ind].velCM[1]=(float)vy/cont;
	H[ind].velCM[2]=(float)vz/cont;
	H[ind].num=cont;
	
}

/*Función que calcula los centros de masa y las velocidades de centro de masa
 * de los halos encontrados*/
void centros(char *archivo){
	FILE *in, *out;
	int i, basura;
	//char archivo[20] = "datos_halos.data";
		
	/*Se abre el archivo con la lista de grupos*/
	if(!(in=fopen(nombregrupos,"r"))){
		printf("Hubo un problema abriendo el archivo %s",nombregrupos);
		exit(1);
	}
	
	/*Se asigna la memoria para los grupos*/
	if(!(G=malloc(NumPart*sizeof(int)))){
		printf("Hubo un problema en la asignación de memoria en los Grupos\n");
		exit(1);
	}
	
	/*Se lee el archivo de lista de grupos, saltándonos el encabezado*/
	fscanf(in,"%d\n",&basura);
	for(i=1;i<=NumPart;i++){
		fscanf(in,"%d\n",&G[i]);
		if(G[i]>max){
			max=G[i];
		}
	}
	fclose(in);	
	printf("The number of halos is %d\n",max);
	
	/*Se asigna la memoria para la tabla de halos*/
	if(!(H=malloc(max*sizeof(struct halo)))){
		printf("Hubo un problema en la asignación de memoria en H\n");
		exit(1);
	}
	
	/*Se procede a hallar centro de masa y velocidad de cada grupo*/
	printf("Calculating centers of mass and velocities\n");
	for(i=1;i<=max;i++){
		datos_halo(i);
	}
	
	/*Se exportan los datos de los halos*/
	printf("Writing Halo data file\n");
	if(!(out=fopen(archivo,"w"))){
		printf("Hubo un problema abriendo el archivo %s",archivo);
		exit(1);
	}
	fprintf(out,"%d\n",max);
	for(i=1;i<=max;i++){
		fprintf(out,"%f %f %f %f %f %f %f\n",H[i].num*header1.mass[1]*pow(10,10),H[i].posCM[0],H[i].posCM[1],H[i].posCM[2],H[i].velCM[0],H[i].velCM[1],H[i].velCM[2]);
	}
	fclose(out);
}

/*Función que realiza el algoritmo de Cloud in cell para
 * calcular la distribución de densidades en la caja*/
void cic(int div){
	float l,x,y,z,fx,fy,fz,wx,wy,wz;
	int k,cx,cy,cz;
	l=header1.BoxSize/div;
	densidades =malloc(sizeof(float)*pow(div,3));
	if(!densidades){
		printf("Hubo un problema con la asignación de memoria\n");
		exit(1);
	}
	printf("The size of the box is %0.1f kPc\nThe lenght of division is %0.1f kPc\n",header1.BoxSize,l);
	for(k=1;k<=NumPart;k++){
		x=P[k].Pos[0];
		y=P[k].Pos[1];
		z=P[k].Pos[2];
				
		cx= (int) x/l;
		cy= (int) y/l;
		cz= (int) z/l;
		
		fx=(x/l-cx)-0.5;
        fy=(y/l-cy)-0.5;
        fz=(z/l-cz)-0.5;
        
        wx=-fabsf(fx)+1;
        wy=-fabsf(fy)+1;
        wz=-fabsf(fz)+1;
        
        add_dens(cx,cy,cz,wx*wy*wz,div);
        add_dens(cx+sign(fx),cy,cz,(1-wx)*wy*wz,div);
        add_dens(cx,cy+sign(fy),cz,wx*(1-wy)*wz,div);
        add_dens(cx,cy,cz+sign(fz),wx*wy*(1-wz),div);
        add_dens(cx+sign(fx),cy+sign(fy),cz,(1-wx)*(1-wy)*wz,div);
        add_dens(cx+sign(fx),cy,cz+sign(fz),(1-wx)*wy*(1-wz),div);
        add_dens(cx,cy+sign(fy),cz+sign(fz),wx*(1-wy)*(1-wz),div);
        add_dens(cx+sign(fx),cy+sign(fy),cz+sign(fz),(1-wx)*(1-wy)*(1-wz),div);    

	}
	/*
	printf("%f\n",pow(l,3));
	for(k=0;k<pow(div,3);k++){
		
		densidades[k]/=pow(l,3);
		printf("%f\n",densidades[k]);
	}*/
	printf("Density field done\n");
	
}

/*Función que recibe el apuntador de densidades, las coordenadas y el peso y asigna el valor establecido*/
void add_dens(int x, int y, int z, float peso, int div){
	x=(x+div)%div;
	y=(y+div)%div;
	z=(z+div)%div;
	densidades[(int)pow(div,2)*x+(int)div*y+z]+=peso*header1.mass[1];
	
}

/*Función Signo*/
int sign(int v)
{
	return v > 0 ? 1 : (v < 0 ? -1 : 0);
}
/*Función que escribe el archivo con las densidades*/
void escribir_densidades(char nom_dens[], int div){
	FILE *Fdens;
	int ind;
	float l=header1.BoxSize/div;
	Fdens=fopen(nom_dens,"w");
	if(!Fdens){
		printf("Hubo problemas abriendo el archivo %s\n",nom_dens);
	}
	
	printf("Writing density file %s\n",nom_dens);
	printf("The file has a single row of data of length %d^3\n",div);
	printf("The format is row-mayor order, being z coordinate the fastest changing index\n");
	printf("The first row is the length of the side of the division in kPc\n");
	
	fprintf(Fdens,"%f\n",l);
	for(ind=0;ind<pow(div,3);ind++){
		fprintf(Fdens,"%f\n",densidades[ind]);
	}
	printf("Density file written\n");
}
