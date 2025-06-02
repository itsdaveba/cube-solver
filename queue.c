#include <stdlib.h>
#include <stdio.h>
#include "queue.h"

static Node *create_node(int depth, long index)
{
    Node *node = malloc(sizeof(Node));
    node->depth = depth;
    node->index = index;
    node->next = NULL;
    return node;
}

Queue *create_queue()
{
    Queue *queue = malloc(sizeof(Queue));
    queue->front = NULL;
    queue->rear = NULL;
    return queue;
}

int is_empty(Queue *queue)
{
    return queue->front == NULL;
}

void enqueue(Queue *queue, int depth, long index)
{
    Node *node = create_node(depth, index);
    if (is_empty(queue))
    {
        queue->front = queue->rear = node;
        return;
    }
    queue->rear->next = node;
    queue->rear = node;
}

Node *dequeue(Queue *queue)
{
    if (is_empty(queue))
        return NULL;
    Node *node = queue->front;
    queue->front = queue->front->next;
    if (queue->front == NULL)
        queue->rear = NULL;
    return node;
}

void free_queue(Queue *queue)
{
    Node *node = queue->front;
    while (node != NULL)
    {
        Node *tmp = node;
        node = node->next;
        free(tmp);
    }
    free(queue);
}

void print_queue(Queue *queue)
{
    Node *node = queue->front;
    while (node != NULL)
    {
        printf("Depth: %d Index: %ld\n", node->depth, node->index);
        node = node->next;
    }
}

int queue_len(Queue *queue)
{
    int len = 0;
    Node *node = queue->front;
    while (node != NULL)
    {
        node = node->next;
        len++;
    }
    return len;
}