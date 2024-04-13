#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <omp.h>

#define MAX_WORD_LENGTH 100

int isWordInList(char *word, char **ppWordList, int numWords)
{
    // check if inputted word is in the ist
    for (int i = 0; i < numWords; i++)
    {
        if (strcmp(word, ppWordList[i]) == 0)
        {
            return 1;
        }
    }
    return 0;
}

int writeListToFile(char **ppList, int length)
{
    char writeFile[32] = {0};

    sprintf(writeFile, "FileOutput.txt");
    FILE *file = fopen(writeFile, "w");
    for (int j = 0; j < length; j++)
    {
        fprintf(file, "%s\n", ppList[j]);
    }
    fclose(file);

    return 0;
}

void writeFileToList(char *fileName, int startIndex, int fileLength, char **ppWordList)
{
    char word[MAX_WORD_LENGTH];

    FILE *file = fopen(fileName, "r"); // pointer to the open file in read, 1 pointer because it is a list of chars

    for (int i = startIndex; i < fileLength; i++)
    {
        fscanf(file, "%s", word);
        ppWordList[i] = strdup(word);
    }

    fclose(file);
}

int readWordFile(char *fileName[], int fileLength[], char ***pppUniqueWordList, int numFiles)
{
    // total length of all lists
    int sumFileLengths = 0;
    for (int i = 0; i < numFiles; i++)
    {
        sumFileLengths += fileLength[i];
    }

    char **ppWordList = (char **)malloc(sumFileLengths * sizeof(char *)); // A 2D list of words being a char ** with each row being a  list of chars to form a string.
    char **ppOutputUniqueWordList = (char **)malloc(sumFileLengths * sizeof(char *));
    int numUniqueWords = 0;

    int startIndex = 0;
    int incFileLength = 0;

    // concatenate all words into one list
    for (int i = 0; i < numFiles; i++)
    {
        // increment where to finish inserting words
        incFileLength += fileLength[i];

        writeFileToList(fileName[i], startIndex, incFileLength, ppWordList);

        // increment where to start inserting words
        startIndex += fileLength[i];
    }

    // iterate through the file to find the unique words
    int i;
#pragma omp parallel for private(i) shared(ppOutputUniqueWordList, ppWordList) schedule(dynamic) num_threads(3)
    for (i = 0; i < sumFileLengths; i++)
    {
        // make word[i] lowercase
        for (int j = 0; ppWordList[i][j]; j++)
        {
            ppWordList[i][j] = tolower(ppWordList[i][j]);
        }

        // check if word is unqiue
        if (!isWordInList(ppWordList[i], ppOutputUniqueWordList, numUniqueWords))
#pragma omp critical
        {
            ppOutputUniqueWordList[numUniqueWords] = strdup(ppWordList[i]);
            numUniqueWords++;
        }
    }

    // reallocate the words into a compact list
    ppOutputUniqueWordList = realloc(ppOutputUniqueWordList, numUniqueWords * sizeof(char *));

    *pppUniqueWordList = ppOutputUniqueWordList;

    return numUniqueWords;
}

unsigned int hash1(const char *str, int len)
{
    unsigned int hash = 0;

    for (int i = 0; i < len; str++, i++)
    {
        hash = (hash + (*str)) << 3;
        str++;
    }

    return hash;
}

unsigned int hash2(const char *str, int len)
{
    unsigned int hash = 0;

    for (int i = 0; i < len; str++, i++)
    {
        hash = (hash << 5) ^ (*str);
        str++;
    }

    return hash;
}

unsigned int hash3(const char *str, int len)
{
    unsigned int hash = 0;

    for (int i = 0; i < len; str++, i++)
    {
        hash ^= (hash << 5) ^ (*str) * (hash + *str);
        str++;
    }

    return hash;
}

unsigned int hash4(const char *str, int len)
{
    unsigned int hash = 0xAAAAAAAA;

    for (int i = 0; i < len; str++, i++)
    {
        hash ^= (hash << 17) ^ ((*str) + (hash & *str));
        str++;
    }

    return hash;
}

void setBit(int *bitArray, int index)
{
    int bitIndex = index / 32;
    int leftShift = index % 32;

    bitArray[bitIndex] |= (1 << leftShift);
}

int getBit(int *bitArray, int index)
{
    int bitIndex = index / 32;
    int leftShift = index % 32;

    return (bitArray[bitIndex] & (1 << leftShift)) != 0;
}

void setbitArray(char **ppUniqueWordList, int *bitArray, int numUniqueWords, int m)
{
    int lenWord;
    int h1, h2, h3, h4;

#pragma omp parallel for private(lenWord, h1, h2, h3, h4) schedule(dynamic) num_threads(3)
    for (int i = 0; i < numUniqueWords; i++)
    {
        // Find the length of the word
        lenWord = 0;
        for (int j = 0; ppUniqueWordList[i][j]; j++)
        {
            lenWord += 1;
        }

        // run hash functions
        h1 = hash1(ppUniqueWordList[i], lenWord) % m;
        h2 = hash2(ppUniqueWordList[i], lenWord) % m;
        h3 = hash3(ppUniqueWordList[i], lenWord) % m;
        h4 = hash4(ppUniqueWordList[i], lenWord) % m;

        // change bits in bit array to 0 indexed by hash functions

        setBit(bitArray, h1);
        setBit(bitArray, h2);
        setBit(bitArray, h3);
        setBit(bitArray, h4);
    }
}

void queryBitArray(int *bitArray, char **ppQueryList, int queryLength, char **insertedList, int insertedListLength, int m)
{
    int lenWord;
    int numFound = 0;
    int numNotFound = 0;
    int falsePositiveCount = 0;
    int bitCountChecker;
    int hashInts[4];

    FILE *file = fopen("query_answer.txt", "w");
#pragma omp parallel for private(lenWord, bitCountChecker, hashInts) reduction(+ : numFound) reduction(+ : numNotFound) reduction(+ : falsePositiveCount) schedule(dynamic) num_threads(3)
    for (int i = 0; i < queryLength; i++)
    {
        lenWord = 0;
        for (int j = 0; ppQueryList[i][j]; j++)
        {
            lenWord += 1;
        }

        // hash each word from query list and check if that position in the bit array is 1. If any positions are 0 it is not in the bitArray.
        hashInts[0] = hash1(ppQueryList[i], lenWord) % m;
        hashInts[1] = hash2(ppQueryList[i], lenWord) % m;
        hashInts[2] = hash3(ppQueryList[i], lenWord) % m;
        hashInts[3] = hash4(ppQueryList[i], lenWord) % m;

        bitCountChecker = 0;

        for (int j = 0; j < 4; j++)
        {
            bitCountChecker += getBit(bitArray, hashInts[j]);
        }

        if (bitCountChecker == 4)
        {
            numFound += 1;
            fprintf(file, "%s 1\n", ppQueryList[i]);
            if (!isWordInList(ppQueryList[i], insertedList, insertedListLength))
            {
                falsePositiveCount += 1;
            }
        }
        else
        {
            numNotFound += 1;
            fprintf(file, "%s 0\n", ppQueryList[i]);
        }
    }

    fclose(file);

    printf("Number of words found: %d\n", numFound);
    printf("Number of words NOT found: %d\n", numNotFound);
}

int main()
{
    struct timespec start, end, startComp, endComp;
    double time_taken;

    clock_gettime(CLOCK_MONOTONIC, &startComp);

    // read files to find unique words
    char *pFiles[4] = {"LITTLE_WOMEN.txt", "MOBY_DICK.txt", "SHAKESPEARE.txt", "MIDDLEMARCH.txt"}; // single pointer as it is a 1D array
    int fileLengths[4] = {195467, 215724, 965465, 319402};

    char **ppUniqueWordList = NULL;
    int numUniqueWords = 0;

    clock_gettime(CLOCK_MONOTONIC, &start);

    numUniqueWords = readWordFile(pFiles, fileLengths, &ppUniqueWordList, 4);

    clock_gettime(CLOCK_MONOTONIC, &end);
    time_taken = (end.tv_sec - start.tv_sec) * 1e9;
    time_taken = (time_taken + (end.tv_nsec - start.tv_nsec)) * 1e-9;
    printf("Total unique words from read files: %d. Process time(s): %lf\n", numUniqueWords, time_taken);

    // Start Bloom Filter Algorithm //

    // initialise bit array
    float falsePositiveRate = 0.05;
    int m = (int)(-(numUniqueWords * log(falsePositiveRate)) / pow(log(2), 2));
    // printf("optimal bit-array length %d\n", m);

    int *BloomFilter = (int *)malloc(m * sizeof(int));

    // int k = (int)((m / numUniqueWords) * log(2));

    // printf("optimal number of hash functions %d\n", k);

    // input unique words into bit array of size n
    clock_gettime(CLOCK_MONOTONIC, &start);

    setbitArray(ppUniqueWordList, BloomFilter, numUniqueWords, m);

    clock_gettime(CLOCK_MONOTONIC, &end);
    time_taken = (end.tv_sec - start.tv_sec) * 1e9;
    time_taken = (time_taken + (end.tv_nsec - start.tv_nsec)) * 1e-9;
    printf("Total process time to set bit array: %lf\n", time_taken);

    int queryLength = 91640;
    char **ppQueryList = (char **)malloc(queryLength * sizeof(char *));

    // read in query file assuming each word is lowercase already
    writeFileToList("query.txt", 0, queryLength, ppQueryList);

    // check if the query words are in the bit array
    clock_gettime(CLOCK_MONOTONIC, &start);

    queryBitArray(BloomFilter, ppQueryList, queryLength, ppUniqueWordList, numUniqueWords, m);

    clock_gettime(CLOCK_MONOTONIC, &end);
    time_taken = (end.tv_sec - start.tv_sec) * 1e9;
    time_taken = (time_taken + (end.tv_nsec - start.tv_nsec)) * 1e-9;
    printf("Total process time to query bit array: %lf\n", time_taken);

    // freeing memory
    for (int i = 0; i < numUniqueWords; i++)
    {
        free(ppUniqueWordList[i]);
    }

    free(ppUniqueWordList);

    free(BloomFilter);

    clock_gettime(CLOCK_MONOTONIC, &endComp);
    time_taken = (endComp.tv_sec - startComp.tv_sec) * 1e9;
    time_taken = (time_taken + (endComp.tv_nsec - startComp.tv_nsec)) * 1e-9;
    printf("Total time taken: %lf\n", time_taken);

    return 0;
}