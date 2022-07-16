import tensorflow as tf
from tensorflow.keras.layers import Dense
from tensorflow.keras.models import Sequential
import numpy as np

real_w = np.random.random((1,5))

real_b = np.random.random(1,)
print(real_w)
print(real_b)
rand_vectors = np.random.random((1000, 5, 1))
z = np.zeros((1,5)).dot(np.zeros((5,1)))
print(f"z shape: {z.shape}")
positive_examples = (real_w[None, :, :].dot(rand_vectors[:500, :, :]) + real_b + np.random.random((1,))).squeeze(axis=-1)
negative_examples = real_w[None, :, :].dot(rand_vectors[500:, :, :]) + real_b - np.random.random((1,))
# first_res_expected = real_w.dot(rand_vectors[349]) + real_b
# first_res = positive_examples[0,0,349,0]
# print(f"equal: {first_res == first_res_expected}")
print(positive_examples.shape)
print(negative_examples.shape)

