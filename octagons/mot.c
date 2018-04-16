#define _GNU_SOURCE
#include <unistd.h>
#include <mqueue.h>
#include <sys/syscall.h>
#include <errno.h>
#include <string.h>
#include <pthread.h>
#include <stdio.h>
#include <semaphore.h>

#define LOOPS	   10000

pthread_t thread[32];

void* worker(void* params);

sem_t sem_start;
sem_t sem_finish;

int main()
{

    sem_init(&sem_start, 0, 0);
    sem_init(&sem_finish, 0, 0);

    for (int k = 0; k < 24; k++)
    {

    	for (int i = 0; i < k; i++)
    	{
		pthread_create(&thread[i], NULL, (void*) worker, NULL);
		pthread_detach(thread[i]);
    	}

    	clock_t begin = clock();

    	for (int j = 0; j < LOOPS; j++)
	{	
    		for(int i=0; i< k; i++) sem_post(&sem_start);

    		for (int i = 0; i < k; i++) sem_wait(&sem_finish);
    	}

    	clock_t end = clock();
    	double delta = ((double)(end - begin) / CLOCKS_PER_SEC) / LOOPS;	
    	printf("threads = %d, mean = %f\n", k, delta);
    }
}    

void* worker(void* params)
{
	while(1) {
	  sem_wait(&sem_start);
	  sem_post(&sem_finish);
	}
}

