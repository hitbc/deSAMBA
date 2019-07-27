#ifndef KTHREAD_H
#define KTHREAD_H

#ifdef __cplusplus
extern "C" {
#endif

void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);
//void * func (void* shared, int step, int tid, void* data);
void kt_pipeline(int n_threads, void *(*func)(void*, int, int, void*), void *shared_data, int n_steps);

#ifdef __cplusplus
}
#endif

#endif
