#ifndef TypeDef_h
#define TypeDef_h

typedef struct yyDetect {
    int channel;
    int ring;
    int sector;
    int energyRaw;
    float energy;
} yyDetect;

typedef struct s3Detect {
    int ring;
    int sector;
    int energyRawRing;
    int energyRawSector;
    float energyRing;
    float energySector;
} s3Detect;

#endif //TypeDef_h