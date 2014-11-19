#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

char nombregrupos[10]="fof.grp";
char directorio[100]="../../outcentro/datos/fof.grp";
int *Ga, *Gc, numPartA, numPartC, numHalosA=0, numHalosC=0, *CorrA, *CorrC;
struct halo{
	int num;
	float posCM[3];
	float velCM[3];
}*Centro,*Act;


void main(int argc, char **argv){
	FILE *ina, *inc;
	int i, basura, j, conteo, mayor, indmay;
	char actual[10]="";
	
	if(argc!=2){
		printf("Debe especificar si se está corriendo desde el directorio mas o menos\n");
		exit(1);
	}
	strcpy(actual,argv[1]);
	
	/*Se abre el archivo con la lista de grupos*/
	if(!(ina=fopen(nombregrupos,"r"))){
		printf("Hubo un problema abriendo el archivo %s",nombregrupos);
		exit(1);
	}
	fscanf(ina,"%d\n",&numPartA);
	
	/*Se asigna la memoria para los grupos*/
	if(!(Ga=malloc(numPartA*sizeof(int)))){
		printf("Hubo un problema en la asignación de memoria en los Grupos\n");
		exit(1);
	}
	
	for(i=0;i<numPartA;i++){
		fscanf(ina,"%d\n",&Ga[i]);
		if(Ga[i]>numHalosA){
			numHalosA=Ga[i];
		}
	}
	fclose(ina);
	
	/*Se abre el archivo con la lista de grupos*/
	if(!(inc=fopen(directorio,"r"))){
		printf("Hubo un problema abriendo el archivo %s",directorio);
		exit(1);
	}
	fscanf(inc,"%d\n",&numPartC);
	
	/*Se asigna la memoria para los grupos*/
	if(!(Gc=malloc(numPartC*sizeof(int)))){
		printf("Hubo un problema en la asignación de memoria en los Grupos\n");
		exit(1);
	}
	
	for(i=0;i<numPartC;i++){
		fscanf(inc,"%d\n",&Gc[i]);
		if(Gc[i]>numHalosC){
			numHalosC=Gc[i];
		}
	}
	fclose(inc);
	/*Se asigna la memoria para los vectores de correspondencia*/
	if(!(CorrA=malloc(numPartA*sizeof(int)))){
		printf("Hubo un problema en la asignación de memoria en los vectores de corresp.\n");
		exit(1);
	}
	
	if(!(CorrC=malloc(numPartC*sizeof(int)))){
		printf("Hubo un problema en la asignación de memoria en los vectores de corresp.\n");
		exit(1);
	}
	
	/*Se llena el vector de correspondencias*/
	for(i=0;i<numHalosC;i++){
		/*Se pone en ceros el vector secundario*/
		for(j=0;j<numHalosA;j++){
			CorrA[j]=0;
		}
		mayor=0;
		indmay=0;
		conteo=0;
		for(j=0;j<numPartC;j++){
			if(Gc[j]==i+1){
				CorrA[Ga[j]-1]++;
				conteo++;
				if(CorrA[Ga[j]-1]>mayor){
					mayor=CorrA[Ga[j]-1];
					indmay=Ga[j]-1;
				}
			}
		}
		CorrC[i]=indmay;
		if(mayor<0.5*conteo){
			printf("El halo %d tiene %d partículas y la mayor correspondencia es con %d con %d partículas\n",i+1,conteo,indmay+1,mayor);
		}
	}
	free(CorrA);
	free(Ga);
	free(Gc);
	
	/*Se importan los datos de los halos*/
	
		
}
