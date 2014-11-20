#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

char nombregrupos[10]="fof.grp";
char directorio[100]="../../outcentro/datos/fof.grp";
char base[40]="halos_";
int *Ga, *Gc, numPartA, numPartC, numHalosA=0, numHalosC=0, *CorrA, *CorrC,na,nc, num;
struct halo{
	float masa;
	float posCM[3];
	float velCM[3];
}*Centro,*Act;

void main(int argc, char **argv){
	FILE *ina, *inc, *inha, *inhc, *out, *indeC, *indeA;
	int i, basura, j, conteo, mayor, indmay;
	char actual[10]="", exportar[50]="centro_a_";
	float deltar, deltav, deltam,rho, l=0.0, *densC, *densA;;
	
	if(argc!=2){
		printf("Debe especificar si se está corriendo desde el directorio mas o menos\n");
		exit(1);
	}
	strcpy(actual,argv[1]);
	
	
	/*Se abre el archivo con la lista de grupos*/
	if(!(ina=fopen(nombregrupos,"r"))){
		printf("Hubo un problema abriendo el archivo %s\n",nombregrupos);
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
		printf("Hubo un problema abriendo el archivo %s\n",directorio);
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
	printf("Realizando correspondencias\n");
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
				num=Ga[j];
				if(num>0){
					CorrA[Ga[j]-1]++;
				}else{
					CorrA[Ga[j]]++;
				}
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
	//free(CorrA);
	//free(Ga);
	//free(Gc);
	printf("Correspondencias listas\n");
	/*Se separa memoria para la información de los halo y se importan sus halos*/
	
	if(!(Centro=malloc(numHalosC*sizeof(struct halo)))){
		printf("Hubo un problema en la asignación de memoria en H\n");
		exit(1);
	}
	if(!(Act=malloc(numHalosA*sizeof(struct halo)))){
		printf("Hubo un problema en la asignación de memoria en H\n");
		exit(1);
	}
	
	if(!(inhc=fopen("../../outcentro/datos/halos_centro","r"))){
		printf("Hubo un problema abriendo el archivo %s\n","halos_centro");
		exit(1);
	}
	printf("Importando información de halos\n");
	fscanf(inhc,"%d\n",&nc);
	if(nc!=numHalosC){
		printf("El número de halos centro no coincide\n");
		exit(1);
	}
	for(i=0;i<numHalosC;i++){
		fscanf(inhc, "%f %f %f %f %f %f %f\n",&Centro[i].masa,&Centro[i].posCM[0],&Centro[i].posCM[1],&Centro[i].posCM[2],&Centro[i].velCM[0],&Centro[i].velCM[1],&Centro[i].velCM[2]);
	}
	fclose(inhc);
	
	strcat(base, actual);
	if(!(inha=fopen(base,"r"))){
		printf("Hubo un problema abriendo el archivo %s\n",base);
		exit(1);
	}
	fscanf(inha,"%d\n",&na);
	if(na!=numHalosA){
		printf("El número de halos centro no coincide\n");
		exit(1);
	}
	for(i=0;i<numHalosA;i++){
		fscanf(inha, "%f %f %f %f %f %f %f\n",&Act[i].masa,&Act[i].posCM[0],&Act[i].posCM[1],&Act[i].posCM[2],&Act[i].velCM[0],&Act[i].velCM[1],&Act[i].velCM[2]);
	}
	fclose(inha);
	printf("Información de Halos importada\n");
	/*Se procede a asignar memoria e importar los archivos de densidades*/
	if(!(densC=malloc(pow(500,3)*sizeof(float)))){
		printf("Hubo un problema en la asignación de memoria en densC\n");
		exit(1);
	}
	if(!(indeC=fopen("../../outcentro/datos/densidades_centro","r"))){
		printf("Hubo un problema abriendo el archivo %s\n","../../outcentro/datos/densidades_centro");
		exit(1);
	}
	fscanf(indeC,"%f\n",&l);
	for(j=0;j<pow(500,3);j++){
		fscanf(indeC,"%f\n", &densC[j]);
	}
	fclose(indeC);
	printf("Se comenzará a exportar el archivo comparativo\n");
	printf("El formato es Masa de referencia, Delta pos, Delta Vel, Delta M\n");
	strcat(exportar,actual);
	if(!(out=fopen(exportar,"r"))){
		printf("Hubo un problema abriendo el archivo %s\n",exportar);
		exit(1);
	}
	for(i=0;i<numHalosC;i++){
		deltar=sqrt(pow(Centro[i].posCM[0]-Act[CorrC[i]].posCM[0],2)+pow(Centro[i].posCM[1]-Act[CorrC[i]].posCM[1],2)+pow(Centro[i].posCM[2]-Act[CorrC[i]].posCM[2],2));
		deltav=sqrt(pow(Centro[i].velCM[0]-Act[CorrC[i]].velCM[0],2)+pow(Centro[i].velCM[1]-Act[CorrC[i]].velCM[1],2)+pow(Centro[i].velCM[2]-Act[CorrC[i]].velCM[2],2));
		deltam=fabsf(Centro[i].masa-Act[i].masa);
		fprintf(out,"%f %f %f %f %f\n",Centro[i].masa,deltar,deltav,deltam,densC[i]);
	}
	fclose(out);
	
}
