#include "read_par_eos.h"

void read_par_eos(char *parfile){
    FILE *fp = fopen(parfile, "r");
    if(fp==NULL){printf("Error opening eos.par\n"); exit(1);}
    const unsigned MAX_LENGTH = 256;
    char buffer[MAX_LENGTH], line2[MAX_LENGTH];

    fgets(buffer, MAX_LENGTH, fp); //skip header
    
    //Read gamma_0
    fgets(buffer, MAX_LENGTH, fp);
    memcpy(line2, buffer, MAX_LENGTH);
    (eos).g0=strtod(line2,NULL);

    fgets(buffer, MAX_LENGTH, fp); //skip line
    fgets(buffer, MAX_LENGTH, fp); //skip header
    //Read gamma_1
    fgets(buffer, MAX_LENGTH, fp);
    memcpy(line2, buffer, MAX_LENGTH);
    (eos).g1=strtod(line2,NULL);

    fgets(buffer, MAX_LENGTH, fp); //skip line
    fgets(buffer, MAX_LENGTH, fp); //skip header
    //Read sigma
    fgets(buffer, MAX_LENGTH, fp);
    memcpy(line2, buffer, MAX_LENGTH);
    (eos).s0=strtod(line2,NULL);

    fgets(buffer, MAX_LENGTH, fp); //skip line
    fgets(buffer, MAX_LENGTH, fp); //skip header
    //Read rho_0
    fgets(buffer, MAX_LENGTH, fp);
    memcpy(line2, buffer, MAX_LENGTH);
    (eos).r0=strtod(line2,NULL);


    fclose(fp);
}




    

    

    
