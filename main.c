#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct sha_256 {
    unsigned int block[64];
    unsigned int K[64];
    unsigned int H[8];
    unsigned int a, b, c, d, e, f, g, h, T1, T2;
    char digest[65];
} Sha256;

void sha256KInit(Sha256 *sha256) {
    int i;
    unsigned int kAux[64] = {
            0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
            0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
            0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
            0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
            0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
            0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
            0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
            0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
    };

    for (i = 0; i < 64; i++)
        (sha256->K)[i] = kAux[i];
}

void sha256HInit(Sha256 *sha256) {
    int i;
    unsigned int hAux[8] = {
            0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
    };

    for (i = 0; i < 8; i++)
        sha256->H[i] = hAux[i];
}

unsigned int ROTR(unsigned int n, unsigned int x) {
    return (x >> n) | (x << (32 - n));
}

unsigned int SHR(int n, unsigned int x) {
    return x >> n;
}

unsigned int Ch(unsigned int x, unsigned int y, unsigned int z) {
    return (x & y) ^ (~x & z);
}

unsigned int Maj(unsigned int x, unsigned int y, unsigned int z) {
    return (x & y) ^ (x & z) ^ (y & z);
}

unsigned int E0(unsigned int x) {
    return ROTR(2, x) ^ ROTR(13, x) ^ ROTR(22, x);
}

unsigned int E1(unsigned int x) {
    return ROTR(6, x) ^ ROTR(11, x) ^ ROTR(25, x);
}

unsigned int sigma0(unsigned int x) {
    return ROTR(7, x) ^ ROTR(18, x) ^ SHR(3, x);
}

unsigned int sigma1(unsigned int x) {
    return ROTR(17, x) ^ ROTR(19, x) ^ SHR(10, x);
}

void printHs(unsigned int *H) {
    int i;
    printf("\n");
    for (i = 0; i < 8; i++) {
        printf("H[%d] = %x\n", i, H[i]);
    }
}

void printBlocks(Sha256 *sha256, int blocksCount) {
    int blockIndex;

    printf("\n\nBlock[%d]\n\n", blocksCount);
    for (blockIndex = 0; blockIndex < 64; blockIndex++)
        printf("block[%d] = %u\n", blockIndex, sha256->block[blockIndex]);
}

void clearBlockLines(unsigned int *block, int start) {
    int index;

    for (index = start; index < 16; index++)
        block[index] = 0;
}

void sha256Encode(char *message, unsigned int msgLength, Sha256 *sha256) {
    unsigned int msgIndex, blockIndex, blockLine, diff, bitAdded, blocksCount;
    char *msgAux;

    bitAdded = blocksCount = 0;

    while ((blocksCount << 6) <= (msgLength + 9)) {
        blockIndex = 0;
        msgAux = ((blocksCount << 6) > msgLength) ? NULL : &(message[blocksCount << 6]);
        blockLine = msgAux ? (msgAux[0] << 8) : 0;
        msgIndex = 1;

        while ((msgIndex < 64) && msgAux && msgAux[msgIndex]) {
            blockLine += msgAux[msgIndex++];
            if (msgIndex % 4 == 0) {
                sha256->block[blockIndex++] = blockLine;
                blockLine = msgAux[msgIndex++];
            }
            blockLine = blockLine << 8;
        }

        if (!bitAdded && msgIndex < 64) {
            blockLine += (1 << 7);
            diff = (((4 - msgLength % 4) << 3) - 1);

            sha256->block[blockIndex++] = blockLine << (diff - 7);
            bitAdded = 1;
        }

        clearBlockLines(sha256->block, blockIndex);

        if (((blocksCount + 1) << 6) > (msgLength + 9)) {
            sha256->block[14] = (msgLength << 3) - ((int) (msgLength << 3));
            sha256->block[15] = (int) (msgLength << 3);
        }

        sha256->a = sha256->H[0];
        sha256->b = sha256->H[1];
        sha256->c = sha256->H[2];
        sha256->d = sha256->H[3];
        sha256->e = sha256->H[4];
        sha256->f = sha256->H[5];
        sha256->g = sha256->H[6];
        sha256->h = sha256->H[7];

        for (blockIndex = 0; blockIndex < 16; blockIndex++) {
            sha256->T1 = sha256->h + E1(sha256->e) + Ch(sha256->e, sha256->f, sha256->g) + sha256->K[blockIndex] +
                         sha256->block[blockIndex];
            sha256->T2 = E0(sha256->a) + Maj(sha256->a, sha256->b, sha256->c);
            sha256->h = sha256->g;
            sha256->g = sha256->f;
            sha256->f = sha256->e;
            sha256->e = sha256->d + sha256->T1;
            sha256->d = sha256->c;
            sha256->c = sha256->b;
            sha256->b = sha256->a;
            sha256->a = sha256->T1 + sha256->T2;
        }

        for (blockIndex = 16; blockIndex < 64; blockIndex++) {
            sha256->block[blockIndex] =
                    sigma1(sha256->block[blockIndex - 2]) + sha256->block[blockIndex - 7] +
                    sigma0(sha256->block[blockIndex - 15]) + sha256->block[blockIndex - 16];

            sha256->T1 = sha256->h + E1(sha256->e) + Ch(sha256->e, sha256->f, sha256->g) + sha256->K[blockIndex] +
                         sha256->block[blockIndex];
            sha256->T2 = E0(sha256->a) + Maj(sha256->a, sha256->b, sha256->c);
            sha256->h = sha256->g;
            sha256->g = sha256->f;
            sha256->f = sha256->e;
            sha256->e = sha256->d + sha256->T1;
            sha256->d = sha256->c;
            sha256->c = sha256->b;
            sha256->b = sha256->a;
            sha256->a = sha256->T1 + sha256->T2;
        }

        sha256->H[0] = (sha256->H[0] + sha256->a) >> 0;
        sha256->H[1] = (sha256->H[1] + sha256->b) >> 0;
        sha256->H[2] = (sha256->H[2] + sha256->c) >> 0;
        sha256->H[3] = (sha256->H[3] + sha256->d) >> 0;
        sha256->H[4] = (sha256->H[4] + sha256->e) >> 0;
        sha256->H[5] = (sha256->H[5] + sha256->f) >> 0;
        sha256->H[6] = (sha256->H[6] + sha256->g) >> 0;
        sha256->H[7] = (sha256->H[7] + sha256->h) >> 0;

        blocksCount++;
    }

    sprintf(sha256->digest, "%08x%08x%08x%08x%08x%08x%08x%08x", sha256->H[0], sha256->H[1], sha256->H[2], sha256->H[3],
            sha256->H[4], sha256->H[5], sha256->H[6], sha256->H[7]);

}

int main() {
    char *message = NULL, fileLine[900];
    Sha256 *sha256;
    FILE *arq;

    arq = fopen("main.c", "r");
    if (arq == NULL) {
        printf("problemas na leitura do arquivo.");
        return 0;
    }

    fgets(fileLine, 900, arq);
    message = (char *) malloc(sizeof(char) * (strlen(fileLine) + 1));
    strcat(message, fileLine);

    while (!feof(arq)) {
        if (fgets(fileLine, 900, arq)) {
            message = (char *) realloc(message, sizeof(char) * (strlen(fileLine) + strlen(message) + 1));
            strcat(message, fileLine);
        }
    }

    sha256 = (Sha256 *) malloc(sizeof(Sha256));

    sha256KInit(sha256);
    sha256HInit(sha256);

    sha256Encode(message, strlen(message), sha256);

    printf("\n\nDigest = %s\n\n\n", sha256->digest);

    free(sha256);

    return 0;
}
