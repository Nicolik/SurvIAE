from enum import IntEnum

from numpy.random import seed
seed(1)  # keras seed fixing import

import tensorflow
tensorflow.random.set_seed(2)  # tensorflow seed fixing

from keras.layers import Input, Dense, concatenate, BatchNormalization, Dropout
from keras.models import Model


class AutoencoderDimension:
    LAYER_1_ENCODER = None
    LAYER_2_ENCODER = None
    CLASSIFIER = None
    CLASSIFIER_HIDDEN = None


class AutoencoderDimensionXS(AutoencoderDimension):
    LAYER_1_ENCODER = 64
    LAYER_2_ENCODER = 32
    CLASSIFIER = 16
    CLASSIFIER_HIDDEN = 8


class AutoencoderDimensionS(AutoencoderDimension):
    LAYER_1_ENCODER = 128
    LAYER_2_ENCODER = 64
    CLASSIFIER = 32
    CLASSIFIER_HIDDEN = 16


class AutoencoderDimensionM(AutoencoderDimension):
    LAYER_1_ENCODER = 256
    LAYER_2_ENCODER = 128
    CLASSIFIER = 64
    CLASSIFIER_HIDDEN = 32


class AutoencoderDimensionL(AutoencoderDimension):
    LAYER_1_ENCODER = 512
    LAYER_2_ENCODER = 256
    CLASSIFIER = 128
    CLASSIFIER_HIDDEN = 64


class AutoencoderDimensionXL(AutoencoderDimension):
    LAYER_1_ENCODER = 1024
    LAYER_2_ENCODER = 512
    CLASSIFIER = 256
    CLASSIFIER_HIDDEN = 128


class AutoencoderDimensionXXL(AutoencoderDimension):
    LAYER_1_ENCODER = 2048
    LAYER_2_ENCODER = 1024
    CLASSIFIER = 512
    CLASSIFIER_HIDDEN = 256


class AutoencoderSize(IntEnum):
    BASELINE = -1
    XS = 0
    S = 1
    M = 2
    L = 3
    XL = 4
    XXL = 5


AE_DIM_DICT = {
    AutoencoderSize.BASELINE: AutoencoderDimension,
    AutoencoderSize.XS: AutoencoderDimensionXS,
    AutoencoderSize.S: AutoencoderDimensionS,
    AutoencoderSize.M: AutoencoderDimensionM,
    AutoencoderSize.L: AutoencoderDimensionL,
    AutoencoderSize.XL: AutoencoderDimensionXL,
    AutoencoderSize.XXL: AutoencoderDimensionXXL,
}


def build_classifier(input_encoded=32, output=1, dimension=AutoencoderDimension, use_hidden_1=False):
    input_data = Input(shape=(input_encoded,))

    if use_hidden_1:
        encoded_clf = BatchNormalization()(input_data)
        classifier = Dense(dimension.CLASSIFIER, activation='relu', name='hidden_1')(encoded_clf)
        classifier = Dropout(0.3)(classifier)
    else:
        classifier = input_data

    classifier_hidden = Dense(dimension.CLASSIFIER_HIDDEN, activation='relu', name='hidden_2')(classifier)
    classifier_hidden = Dropout(0.3)(classifier_hidden)
    classifier_hidden = BatchNormalization()(classifier_hidden)

    classifier_output = Dense(output, activation='sigmoid', name='clf_layer')(classifier_hidden)

    classifier_model = Model(input_data, classifier_output)
    return classifier_model


def build_classifier_clinical(input_clinical=8, output=1, input_encoded=32, dimension=AutoencoderDimension):
    input_data = Input(shape=(input_encoded,))
    clinical_data = Input(shape=(input_clinical,))

    encoded_clf = BatchNormalization()(input_data)
    classifier = Dense(dimension.CLASSIFIER, activation='relu')(encoded_clf)
    classifier = Dropout(0.3)(classifier)

    merged = concatenate([classifier, clinical_data])
    merged_clf = BatchNormalization()(merged)
    classifier_hidden = Dense(dimension.CLASSIFIER_HIDDEN, activation='relu')(merged_clf)

    classifier_hidden = Dropout(0.3)(classifier_hidden)
    classifier_hidden = BatchNormalization()(classifier_hidden)
    classifier_output = Dense(output, activation='sigmoid', name='clf_layer')(classifier_hidden)

    classifier_model = Model([input_data, clinical_data], classifier_output)

    return classifier_model


def build_autoencoder(dimension=AutoencoderDimension, input_rna=15945):
    # Input layer
    input_data = Input(shape=(input_rna,))

    # Encoder layers
    encoded = Dense(dimension.LAYER_1_ENCODER, activation='relu', name='layer_1_encoder')(input_data)
    encoded = Dense(dimension.LAYER_2_ENCODER, activation='relu', name='layer_2_encoder')(encoded)

    # Decoder layers
    decoded = Dense(dimension.LAYER_1_ENCODER, activation='relu', name='layer_1_decoder')(encoded)
    decoded = Dense(input_rna, activation='relu', name='layer_2_decoder')(decoded)

    # Autoencoder model
    autoencoder = Model(input_data, decoded)

    # Encoder model
    encoder = Model(input_data, encoded)

    # Decoder model
    decoder = Model(encoded, decoded)

    return autoencoder, encoder, decoder


def extract_encoder(autoencoder):
    return Model(
        inputs=autoencoder.input,
        outputs=autoencoder.get_layer('layer_2_encoder').output
    )


def extract_hidden(classifier):
    return Model(
        inputs=classifier.input,
        outputs=classifier.get_layer('hidden_2').output
    )


def extract_encoder_hidden(encoder_classifier):
    clf = encoder_classifier.layers[-1]
    hidden = extract_hidden(clf)
    encoded = encoder_classifier.layers[-2].output
    hidden_out = hidden(encoded)
    return Model(
        inputs=encoder_classifier.input,
        outputs=hidden_out
    )


def extract_decoder(autoencoder):
    return Model(
        inputs=autoencoder.get_layer('layer_2_encoder').output,
        outputs=autoencoder.get_layer('layer_2_decoder').output
    )


def combine_models(encoder, classifier):
    for layer in encoder.layers:
        layer.trainable = False

    encoded = encoder.output
    clf_output = classifier(encoded)

    return Model(
        inputs=encoder.input,
        outputs=clf_output
    )
