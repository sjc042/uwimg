from uwimg import *

def softmax_model(inputs, outputs):
    l = [make_layer(inputs, outputs, SOFTMAX)]
    # l = [   make_layer(inputs, outputs, SOFTMAX),
            # make_layer(64, 32, SOFTMAX),
            # make_layer(32, outputs, SOFTMAX)]
    return make_model(l)

def neural_net(inputs, outputs):
    print(inputs)
    l = [   make_layer(inputs, 32, RELU),
            make_layer(32, outputs, SOFTMAX)]
    # l = [   make_layer(inputs, 64, RELU),
    #     make_layer(64, 32, SOFTMAX),
    #     make_layer(32, outputs, SOFTMAX)]
    return make_model(l)

print("loading data...")
# train = load_classification_data(b"mnist.train", b"mnist.labels", 1)
# test  = load_classification_data(b"mnist.test",  b"mnist.labels", 1)
train = load_classification_data(b"cifar.train", b"cifar/labels", 1)
test  = load_classification_data(b"cifar.test",  b"cifar/labels", 1)
print("done")
print


batch = 128
iters = 3000
rate = 0.001      # originally 0.01
momentum = .9   # originally 0.9
decay = 0.1      # originally 0

nn = neural_net(train.X.cols, train.y.cols)
sm = softmax_model(train.X.cols, train.y.cols)
model_str = 'Neural Net'
model = nn if (model_str == 'Neural Net') else sm
print(f"training model {model_str}...")
train_model(model, train, batch, iters, rate, momentum, decay)
print("done")
print

print("evaluating model...")
print("training accuracy:", accuracy_model(model, train))
print("test accuracy:    ", accuracy_model(model, test))


