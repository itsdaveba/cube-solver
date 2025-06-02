#ifndef QUEUE_H
#define QUEUE_H

typedef struct Node
{
    int depth;
    long index;
    struct Node *next;
} Node;

typedef struct Queue
{
    Node *front;
    Node *rear;
} Queue;

Queue *create_queue();

int is_empty(Queue *queue);

void enqueue(Queue *queue, int depth, long index);

Node *dequeue(Queue *queue);

void free_queue(Queue *queue);

void print_queue(Queue *queue);

int queue_len(Queue *queue);

#endif
