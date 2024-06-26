#include <stdio.h>
#include <string.h>

int main() {
    char dataset[] = "GAGTCGAGAACAACTGCTCGGACCGCATTTCTCTTCATGCCCGGACAATTAACCCAGAATACCCTAACCCCGATGAAGCAGCTAAAGATACAAGGGGCAGTGTGGTATAGATCGTCATCGCCCCTTGAACAAGTAGTGTCCGGAGGGACGTAGTGCCTGCCCGCTGTGTTTACGGTCTTAGTATAGTCTATGACGCGAGTGCCAAGAGCGTCATTGTAGCGCGGTAGAGCCAGCTGGACAGGCTTCTGATAGTACTGGATTACGCTCAGGTGGAATAGTCTTAGACGATATCGGCCCTATTGCCGAATTGACGAGGACGCTAGACTACTGCCACGCGAATAGATCGCTAAGCGGGTAAAGCTAGCAGGGGTTATCGCGACGTTTTGTGAAAGGCGACGATTATCTGACGGGTACTTCCGGCCTTGTAGTTAGACCTGACGCCGGCAGAAAGAAGTATGAACAGTCTTGCCCAGGATACACGCTTTACTTACGCACGCGGCCTATCGACCTGGTGTACTTTCGCCGCTCGCCGGGAACGGGCAACATTAACGTCATCAAAGTTTCCCTCACGGTATCCCCCATTGACACATGAATAACACCGCAGCATCTGACAGGAACAGCTTGTACTGTATTAGGGTGATTAGATTTACCAGCTGGCACGACTCAATCCTTCATCACATTAACCTGTCTGACCTTGGCGACATTAGGGCAGATAAAATGGTCGAATCAGAAAGTAGGCAAGGTTTAGCAAATAACGCACGGTTGGTAGCACGCAAAATTAGTGCAGATATACATGAGCAGGCACATATTAATTCTGTTTGTCGCTTATGCACTTCTTTGTAAGGTCACACAAAGTCAACTTGCCGGCAATGGAGACACCGTTGGGGCATGTTAGAGGGCCGGACGGT";
    
    int A = 0, C = 0, G = 0, T = 0;
    int i;

    for (i = 0; i < strlen(dataset); i++) {
        switch (dataset[i]) {
            case 'A':
                A++;
                break;
            case 'C':
                C++;
                break;
            case 'G':
                G++;
                break;
            case 'T':
                T++;
                break;
        }
    }

    printf("A: %d\n", A);
    printf("C: %d\n", C);
    printf("G: %d\n", G);
    printf("T: %d\n", T);

    return 0;
}
