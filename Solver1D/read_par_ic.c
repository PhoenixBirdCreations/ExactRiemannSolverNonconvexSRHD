#include "read_par_ic.h"

char *strip_copy(char const *s)
{
    char *buf = malloc(1 + strlen(s));
    if (buf)
    {
        char *p = buf;
        char const *q;
        int n;
        for (q = s; *q; q += n + strspn(q+n, "\n"))
        {
            n = strcspn(q, "\n");
            strncpy(p, q, n);
            p += n;
        }
        *p++ = '\0';
        buf = realloc(buf, p - buf);
    }
    return buf;
}

void read_par_ic(char *parfile, double *ciInfo, double *goalInfo, char **exact_solution_file, char **wave_curves_file){
    FILE *fp = fopen(parfile, "r");
    if(fp==NULL){printf("Error opening %s\n",parfile); exit(1);}
    const unsigned MAX_LENGTH = 256;
    char buffer[MAX_LENGTH], line2[MAX_LENGTH];

    fgets(buffer, MAX_LENGTH, fp); //skip header
    
    //Initial conditions
    for(int k=0; k<6; k++){
        fgets(buffer, MAX_LENGTH, fp);
        memcpy(line2, buffer, MAX_LENGTH);
        ciInfo[k]=strtod(line2,NULL);
    }
    

    //goals for pressure and density
    fgets(buffer, MAX_LENGTH, fp); //skip line
    fgets(buffer, MAX_LENGTH, fp); //skip header
    for(int k=0; k<4; k++){
        fgets(buffer, MAX_LENGTH, fp);
        memcpy(line2, buffer, MAX_LENGTH);
        goalInfo[k]=strtod(line2,NULL);
    }

    //final time
    fgets(buffer, MAX_LENGTH, fp); //skip line
    fgets(buffer, MAX_LENGTH, fp); //skip header    
    fgets(buffer, MAX_LENGTH, fp);
    memcpy(line2, buffer, MAX_LENGTH);
    ciInfo[6]=strtod(line2,NULL);

    //output file name
    fgets(buffer, MAX_LENGTH, fp); //skip line
    fgets(buffer, MAX_LENGTH, fp); //skip header
    fgets(buffer, MAX_LENGTH, fp);
    memcpy(line2, buffer, MAX_LENGTH); 
    (*exact_solution_file)=strip_copy(line2);

    //output file name
    fgets(buffer, MAX_LENGTH, fp); //skip line
    fgets(buffer, MAX_LENGTH, fp); //skip header
    fgets(buffer, MAX_LENGTH, fp);
    memcpy(line2, buffer, MAX_LENGTH); 
    (*wave_curves_file)=strip_copy(line2);
    
    fclose(fp);
}


