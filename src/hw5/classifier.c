#include <math.h>
#include <stdlib.h>
#include "image.h"
#include "matrix.h"

// Run an activation function on each element in a matrix,
// modifies the matrix in place
// matrix m: Input to activation function
// ACTIVATION a: function to run
void activate_matrix(matrix m, ACTIVATION a)
{
    int i, j;
    for(i = 0; i < m.rows; ++i){
        double sum = 0;
        for(j = 0; j < m.cols; ++j){
            double x = m.data[i][j];
            if(a == LOGISTIC){
                // TODO
                m.data[i][j] = 1 / (1 + exp(-1 * x));
            } else if (a == RELU){
                // TODO
                m.data[i][j] = (x > 0) ? x : 0;
            } else if (a == LRELU){
                // TODO
                m.data[i][j] = (x > 0) ? x : 0.1 * x;
            } else if (a == SOFTMAX){
                // TODO
                m.data[i][j] = exp(x);
            }
            sum += m.data[i][j];
        }
        if (a == SOFTMAX) {
            // TODO: have to normalize by sum if we are using SOFTMAX
            for(j = 0; j < m.cols; ++j) {
                m.data[i][j] = m.data[i][j] / sum;
            }
        }
    }
}

// Calculates the gradient of an activation function and multiplies it into
// the delta for a layer
// matrix m: an activated layer output
// ACTIVATION a: activation function for a layer
// matrix d: delta before activation gradient
void gradient_matrix(matrix m, ACTIVATION a, matrix d)
{
    int i, j;
    for(i = 0; i < m.rows; ++i){
        for(j = 0; j < m.cols; ++j){
            double x = m.data[i][j];
            // TODO: multiply the correct element of d by the gradient
            double d_act;
            if(a == LOGISTIC){
                // TODO
                d_act = x * (1 - x);
            } else if (a == RELU){
                // TODO
                d_act = (x > 0) ? 1 : 0;
            } else if (a == LRELU){
                // TODO
                d_act = (x > 0) ? 1 : 0.1;
            } else if (a == SOFTMAX){
                // TODO
                d_act = 1;
            }
            d.data[i][j] = d_act * d.data[i][j];
        }
    }
}

// Forward propagate information through a layer
// layer *l: pointer to the layer
// matrix in: input to layer
// returns: matrix that is output of the layer
matrix forward_layer(layer *l, matrix in)
{

    l->in = in;  // Save the input for backpropagation


    // TODO: fix this! multiply input by weights and apply activation function.
    matrix out = matrix_mult_matrix(l->in, l->w);
    activate_matrix(out, l->activation);
    free_matrix(l->out);// free the old output
    l->out = out;       // Save the current output for gradient calculation
    return out;
}

// Backward propagate derivatives through a layer
// layer *l: pointer to the layer
// matrix delta: partial derivative of loss w.r.t. output of layer
// returns: matrix, partial derivative of loss w.r.t. input to layer
matrix backward_layer(layer *l, matrix delta)
{
    // 1.4.1
    // delta is dL/dy
    // TODO: modify it in place to be dL/d(xw)
    // dL/d(xw) = dL/dy * dy/d(xw) = dL/dy * df(xw)/d(xw) = dL/dy * f'(xw) where f(xw) is activation function
    gradient_matrix(l->out, l->activation, delta);

    // 1.4.2
    // TODO: then calculate dL/dw and save it in l->dw
    // dL/dw = dL/d(xw) * d(xw)/dw = dL/d(xw) * x ---> x_transpose * dL/d(xw)
    free_matrix(l->dw);
    matrix dw = matrix_mult_matrix(transpose_matrix(l->in), delta);
    l->dw = dw;

    // 1.4.3
    // TODO: finally, calculate dL/dx and return it.
    // dL/dx = dL/d(xw) * d(xw)/dx = dL/d(xw) * w ---> dL/d(xw) * w_transpose
    matrix dx = matrix_mult_matrix(delta, transpose_matrix(l->w));

    return dx;
}

// Update the weights at layer l
// layer *l: pointer to the layer
// double rate: learning rate
// double momentum: amount of momentum to use
// double decay: value for weight decay
void update_layer(layer *l, double rate, double momentum, double decay)
{
    // TODO:
    // Calculate Δw_t = dL/dw_t - λw_t + mΔw_{t-1}
    // save it to l->v
    matrix v = make_matrix(l->v.rows, l->v.cols);
    int i, j;
    for(i = 0; i < v.rows; ++i) {
        for(j = 0; j < v.cols; ++j) {
            v.data[i][j] = l->dw.data[i][j] - decay * l->w.data[i][j] + momentum * l->v.data[i][j];
        }
    }
    free_matrix(l->v);
    l->v = v;
    // Update l->w
    // w_{t+1} = w_t + ηΔw_t
    matrix w = make_matrix(l->w.rows, l->w.cols);
    for(i = 0; i < w.rows; ++i) {
        for(j = 0; j < w.cols; ++j) {
            w.data[i][j] = l->w.data[i][j] + rate * l->v.data[i][j];
        }
    }
    free_matrix(l->w);
    l->w = w;
    // Remember to free any intermediate results to avoid memory leaks
}

// Make a new layer for our model
// int input: number of inputs to the layer
// int output: number of outputs from the layer
// ACTIVATION activation: the activation function to use
layer make_layer(int input, int output, ACTIVATION activation)
{
    layer l;
    l.in  = make_matrix(1,1);
    l.out = make_matrix(1,1);
    l.w   = random_matrix(input, output, sqrt(2./input));
    l.v   = make_matrix(input, output);
    l.dw  = make_matrix(input, output);
    l.activation = activation;
    return l;
}

// Run a model on input X
// model m: model to run
// matrix X: input to model
// returns: result matrix
matrix forward_model(model m, matrix X)
{
    int i;
    for(i = 0; i < m.n; ++i){
        X = forward_layer(m.layers + i, X);
    }
    return X;
}

// Run a model backward given gradient dL
// model m: model to run
// matrix dL: partial derivative of loss w.r.t. model output dL/dy
void backward_model(model m, matrix dL)
{
    matrix d = copy_matrix(dL);
    int i;
    for(i = m.n-1; i >= 0; --i){
        matrix prev = backward_layer(m.layers + i, d);
        free_matrix(d);
        d = prev;
    }
    free_matrix(d);
}

// Update the model weights
// model m: model to update
// double rate: learning rate
// double momentum: amount of momentum to use
// double decay: value for weight decay
void update_model(model m, double rate, double momentum, double decay)
{
    int i;
    for(i = 0; i < m.n; ++i){
        update_layer(m.layers + i, rate, momentum, decay);
    }
}

// Find the index of the maximum element in an array
// double *a: array
// int n: size of a, |a|
// returns: index of maximum element
int max_index(double *a, int n)
{
    if(n <= 0) return -1;
    int i;
    int max_i = 0;
    double max = a[0];
    for (i = 1; i < n; ++i) {
        if (a[i] > max){
            max = a[i];
            max_i = i;
        }
    }
    return max_i;
}

// Calculate the accuracy of a model on some data d
// model m: model to run
// data d: data to run on
// returns: accuracy, number correct / total
double accuracy_model(model m, data d)
{
    matrix p = forward_model(m, d.X);
    int i;
    int correct = 0;
    for(i = 0; i < d.y.rows; ++i){
        if(max_index(d.y.data[i], d.y.cols) == max_index(p.data[i], p.cols)) ++correct;
    }
    return (double)correct / d.y.rows;
}

// Calculate the cross-entropy loss for a set of predictions
// matrix y: the correct values
// matrix p: the predictions
// returns: average cross-entropy loss over data points, 1/n Σ(-ylog(p))
double cross_entropy_loss(matrix y, matrix p)
{
    int i, j;
    double sum = 0;
    for(i = 0; i < y.rows; ++i){
        for(j = 0; j < y.cols; ++j){
            sum += -y.data[i][j]*log(p.data[i][j]);
        }
    }
    return sum/y.rows;
}


// Train a model on a dataset using SGD
// model m: model to train
// data d: dataset to train on
// int batch: batch size for SGD
// int iters: number of iterations of SGD to run (i.e. how many batches)
// double rate: learning rate
// double momentum: momentum
// double decay: weight decay
void train_model(model m, data d, int batch, int iters, double rate, double momentum, double decay)
{
    int e;
    for(e = 0; e < iters; ++e){
        data b = random_batch(d, batch);
        matrix p = forward_model(m, b.X);
        fprintf(stderr, "%06d: Loss: %f\n", e, cross_entropy_loss(b.y, p));
        matrix dL = axpy_matrix(-1, p, b.y); // partial derivative of loss dL/dy
        backward_model(m, dL);
        update_model(m, rate/batch, momentum, decay);
        free_matrix(dL);
        free_data(b);
    }
}


// Questions 
//
// 5.2.2.1 Why might we be interested in both training accuracy and testing accuracy? What do these two numbers tell us about our current model?
// TODO
// While training accuracy tells us how well the model is fitting the training data, testing accuracy tells us how well the model generalizes by predicting on data it has not seen before.

// 5.2.2.2 Try varying the model parameter for learning rate to different powers of 10 (i.e. 10^1, 10^0, 10^-1, 10^-2, 10^-3) and training the model. What patterns do you see and how does the choice of learning rate affect both the loss during training and the final model accuracy?
// TODO
// A larger learning rate may speed up the gradient descent but may also cause more fluctuations to the loss. In this case, learning rates in 0.1 yield better final test accuracy ~92% while having loss during training around 0 to 1.

// 5.2.2.3 Try varying the parameter for weight decay to different powers of 10: (10^0, 10^-1, 10^-2, 10^-3, 10^-4, 10^-5). How does weight decay affect the final model training and test accuracy?
// TODO
// Changing the decay from 1 to 0 increases both the training and testing accuracies, but probably since parameters are small, any further decrement does not change both accuracies much.

// 5.2.3.1 Currently the model uses a logistic activation for the first layer. Try using a the different activation functions we programmed. How well do they perform? What's best?
// TODO
// Softmax has the best train/test accuracies, overall LRELU(0.95) ~ RELU(0.95) > SOFTMAX(0.94) > LOGISTIC(0.098)

// 5.2.3.2 Using the same activation, find the best (power of 10) learning rate for your model. What is the training accuracy and testing accuracy?
// TODO
// With learning rate = 0.1, training and test accuracy are both ~0.96

// 5.2.3.3 Right now the regularization parameter `decay` is set to 0. Try adding some decay to your model. What happens, does it help? Why or why not may this be?
// TODO
// Changing the decay from 1 to 0 increases both the training and testing accuracies, but since parameters are small, any further decrement causes accuracies to decrease slightly.

// 5.2.3.4 Modify your model so it has 3 layers instead of two. The layers should be `inputs -> 64`, `64 -> 32`, and `32 -> outputs`. Also modify your model to train for 3000 iterations instead of 1000. Look at the training and testing error for different values of decay (powers of 10, 10^-4 -> 10^0). Which is best? Why?
// TODO
// By making the learning rate 0.0001 and decay 0, the neural network reaches 0.30 test accuracy, bigger values of learning yield -nan values in loss.
// I expect the accuracy would improve with a more complex model. Learning rate could be relatively small since we run 3000 iterations. Changes in Decay might not be really significant.

// 5.3.2.1 How well does your network perform on the CIFAR dataset?
// TODO
// The 2-layered neural network reaches about 0.35 test accuracy on the CIFAR dataset, with learning rate 0.001 and decay of 0.
// The execution is killed while running 3-layered network probably due to lack of memory in Virtual Machine.



