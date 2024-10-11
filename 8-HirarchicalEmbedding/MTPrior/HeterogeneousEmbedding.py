import time
import numpy as np
import tensorflow as tf
from jhyexp import my_KNN, my_Kmeans

from models import GAT, HeteGAT, HeteGAT_multi
from utils import process

from CLassificationAnalysis import RClassification, SVMClassification, RFClassification, LRClassification, AdaBoostClassification, NaiveBayesClassification
import os

os.environ["CUDA_VISIBLE_DEVICES"] = "1,2,3"

config = tf.ConfigProto()
config.gpu_options.allow_growth = True

###CHANGE NEEDED
dataset = 'MT'
featype = 'fea'
checkpt_file = '/home/fatemeh/check_point/{}{}_allMP_multi_{}_.ckpt'.format(dataset, dataset, featype)
print('model: {}'.format(checkpt_file))
# training params
batch_size = 1
nb_epochs = 200
patience = 100
lr = 0.005  # learning rate
l2_coef = 0.001  # weight decay
# numbers of hidden units per each attention head in each layer
hid_units = [16] #8  16>256  32>512 ###16
n_heads = [8, 1]  # additional entry for the output layer
residual = False
nonlinearity = tf.nn.elu
model = HeteGAT_multi


print('Dataset: ' + dataset)
print('----- Opt. hyperparams -----')
print('lr: ' + str(lr))
print('l2_coef: ' + str(l2_coef))
print('----- Archi. hyperparams -----')
print('nb. layers: ' + str(len(hid_units)))
print('nb. units per layer: ' + str(hid_units))
print('nb. attention heads: ' + str(n_heads))
print('residual: ' + str(residual))
print('nonlinearity: ' + str(nonlinearity))
print('model: ' + str(model))

# jhy data
import scipy.io as sio
import scipy.sparse as sp



def sample_mask(idx, l):
    """Create mask."""
    mask = np.zeros(l)
    mask[idx] = 1
    return np.array(mask, dtype=np.bool)


import numpy as np
def load_data_dblp():
    in_gmig = np.loadtxt('/home/fatemeh/MTPrior/data/HCC/Import_Genetest/Ready_gmig.txt')
    in_gmicmig = np.loadtxt('/home/fatemeh/MTPrior/data/HCC/Import_Genetest/Ready_gmicmig.txt')
    in_glg = np.loadtxt('/home/fatemeh/MTPrior/data/HCC/Import_Genetest/Ready_glncg.txt')
    in_gmilmig = np.loadtxt('/home/fatemeh/MTPrior/data/HCC/Import_Genetest/Ready_gmilncmig.txt')

    in_gmig2 = np.loadtxt('/home/fatemeh/MTPrior/data/HCC/Import_Genetest/Ready_gmig2_Test.txt')
    in_gmicmig2 = np.loadtxt('/home/fatemeh/MTPrior/data/HCC/Import_Genetest/Ready_gmicmig2_Test.txt')
    in_glg2 = np.loadtxt('/home/fatemeh/MTPrior/data/HCC/Import_Genetest/Ready_glncg2_Test.txt')
    in_gmilmig2 = np.loadtxt('/home/fatemeh/MTPrior/data/HCC/Import_Genetest/Ready_gmilncmig2_Test.txt')

    in_gg1 = np.loadtxt('/home/fatemeh/MTPrior/data/HCC/Import_Genetest/Ready_gg1.txt')
    in_gg2 = np.loadtxt('/home/fatemeh/MTPrior/data/HCC/Import_Genetest/Ready_gg2.txt')

    in_feature = np.loadtxt('/home/fatemeh/MTPrior/data/HCC/Import_Genetest/feature_Exp.txt').astype(float)

    in_label = np.loadtxt('/home/fatemeh/MTPrior/data/HCC/Import_Genetest/Ready_label.txt')
    in_train = np.loadtxt('/home/fatemeh/MTPrior/data/HCC/Import_Genetest/Ready_train.idx2.txt')
    in_val = np.loadtxt('/home/fatemeh/MTPrior/data/HCC/Import_Genetest/Ready_val.idx2.txt')
    in_test = np.loadtxt('/home/fatemeh/MTPrior/data/HCC/Import_Genetest/Ready_test.idx2.txt')





    truelabels, truefeatures = in_label, in_feature.astype(float)
    N = truefeatures.shape[0]
    #rownetworks = [in_gmig - np.eye(N), in_gmicmig - np.eye(N) , in_glg - np.eye(N),
             #in_gmilmig - np.eye(N), in_gg1 - np.eye(N)]

    rownetworks = [in_gmig - np.eye(N), in_gmicmig - np.eye(N), in_glg - np.eye(N), in_gmilmig - np.eye(N),
                    in_gg1 - np.eye(N), in_gmig2 - np.eye(N), in_gmicmig2 - np.eye(N), in_glg2 - np.eye(N), in_gmilmig2 - np.eye(N)]

    y = truelabels

    in_train = np.expand_dims(in_train, axis=0)
    in_val= np.expand_dims(in_val, axis=0)
    in_test = np.expand_dims(in_test, axis=0)


    train_idx = in_train.astype(int)
    val_idx = in_val.astype(int)
    test_idx = in_test.astype(int)

    train_mask = sample_mask(train_idx, y.shape[0])
    val_mask = sample_mask(val_idx, y.shape[0])
    test_mask = sample_mask(test_idx, y.shape[0])

    y_train = np.zeros(y.shape)
    y_val = np.zeros(y.shape)
    y_test = np.zeros(y.shape)

    y_train[train_mask, :] = y[train_mask, :]
    y_val[val_mask, :] = y[val_mask, :]
    y_test[test_mask, :] = y[test_mask, :]


    print('y_train:{}, y_val:{}, y_test:{}, train_idx:{}, val_idx:{}, test_idx:{}'.format(y_train.shape,
                                                                                          y_val.shape,
                                                                                          y_test.shape,
                                                                                          train_idx.shape,
                                                                                          val_idx.shape,
                                                                                          test_idx.shape))
    truefeatures_list = [truefeatures, truefeatures, truefeatures]
    return rownetworks, truefeatures_list, y_train, y_val, y_test, train_mask, val_mask, test_mask , truelabels


adj_list, fea_list, y_train, y_val, y_test, train_mask, val_mask, test_mask, labels = load_data_dblp()
# use adj_list as fea_list, have a try~
if featype == 'adj':
    fea_list = adj_list


import scipy.sparse as sp



nb_nodes = fea_list[0].shape[0]
ft_size = fea_list[0].shape[1]
nb_classes = y_train.shape[1]

print(nb_nodes)
print(ft_size)
print(nb_classes)
# adj = adj.todense()

# features = features[np.newaxis]  # [1, nb_node, ft_size]
fea_list = [fea[np.newaxis] for fea in fea_list]
adj_list = [adj[np.newaxis] for adj in adj_list]
y_train = y_train[np.newaxis]
y_val = y_val[np.newaxis]
y_test = y_test[np.newaxis]
train_mask = train_mask[np.newaxis]
val_mask = val_mask[np.newaxis]
test_mask = test_mask[np.newaxis]

biases_list = [process.adj_to_bias(adj, [nb_nodes], nhood=1) for adj in adj_list]

####add
print('build graph...')
with tf.Graph().as_default():
    with tf.name_scope('input'):
        ftr_in_list = [tf.placeholder(dtype=tf.float32,
                                      shape=(batch_size, nb_nodes, ft_size),
                                      name='ftr_in_{}'.format(i))
                       for i in range(len(fea_list))]
        bias_in_list = [tf.placeholder(dtype=tf.float32,
                                       shape=(batch_size, nb_nodes, nb_nodes),
                                       name='bias_in_{}'.format(i))
                        for i in range(len(biases_list))]
        lbl_in = tf.placeholder(dtype=tf.int32, shape=(
            batch_size, nb_nodes, nb_classes), name='lbl_in')
        msk_in = tf.placeholder(dtype=tf.int32, shape=(batch_size, nb_nodes),
                                name='msk_in')
        attn_drop = tf.placeholder(dtype=tf.float32, shape=(), name='attn_drop')
        ffd_drop = tf.placeholder(dtype=tf.float32, shape=(), name='ffd_drop')
        is_train = tf.placeholder(dtype=tf.bool, shape=(), name='is_train')
    # forward
    logits, final_embedding, att_val = model.inference(ftr_in_list, nb_classes, nb_nodes, is_train,
                                                       attn_drop, ffd_drop,
                                                       bias_mat_list=bias_in_list,
                                                       hid_units=hid_units, n_heads=n_heads,
                                                       residual=residual, activation=nonlinearity)



    # cal masked_loss
    log_resh = tf.reshape(logits, [-1, nb_classes])
    lab_resh = tf.reshape(lbl_in, [-1, nb_classes])
    msk_resh = tf.reshape(msk_in, [-1])
    loss = model.masked_softmax_cross_entropy(log_resh, lab_resh, msk_resh)
    accuracy = model.masked_accuracy(log_resh, lab_resh, msk_resh)
    # optimzie
    train_op = model.training(loss, lr, l2_coef)

    saver = tf.train.Saver()

    init_op = tf.group(tf.global_variables_initializer(),
                       tf.local_variables_initializer())

    vlss_mn = np.inf
    vacc_mx = 0.0
    curr_step = 0

    with tf.Session(config=config) as sess:
        sess.run(init_op)

        train_loss_avg = 0
        train_acc_avg = 0
        val_loss_avg = 0
        val_acc_avg = 0

        for epoch in range(nb_epochs):
            tr_step = 0

            tr_size = fea_list[0].shape[0]
            # ================   training    ============
            while tr_step * batch_size < tr_size:

                fd1 = {i: d[tr_step * batch_size:(tr_step + 1) * batch_size]
                       for i, d in zip(ftr_in_list, fea_list)}
                fd2 = {i: d[tr_step * batch_size:(tr_step + 1) * batch_size]
                       for i, d in zip(bias_in_list, biases_list)}
                fd3 = {lbl_in: y_train[tr_step * batch_size:(tr_step + 1) * batch_size],
                       msk_in: train_mask[tr_step * batch_size:(tr_step + 1) * batch_size],
                       is_train: True,
                       attn_drop: 0.6,
                       ffd_drop: 0.6}
                fd = fd1
                fd.update(fd2)
                fd.update(fd3)
                _, loss_value_tr, acc_tr, att_val_train = sess.run([train_op, loss, accuracy, att_val],
                                                                   feed_dict=fd)
                train_loss_avg += loss_value_tr
                train_acc_avg += acc_tr
                tr_step += 1

            vl_step = 0
            vl_size = fea_list[0].shape[0]
            # =============   val       =================
            while vl_step * batch_size < vl_size:
                # fd1 = {ftr_in: features[vl_step * batch_size:(vl_step + 1) * batch_size]}
                fd1 = {i: d[vl_step * batch_size:(vl_step + 1) * batch_size]
                       for i, d in zip(ftr_in_list, fea_list)}
                fd2 = {i: d[vl_step * batch_size:(vl_step + 1) * batch_size]
                       for i, d in zip(bias_in_list, biases_list)}
                fd3 = {lbl_in: y_val[vl_step * batch_size:(vl_step + 1) * batch_size],
                       msk_in: val_mask[vl_step * batch_size:(vl_step + 1) * batch_size],
                       is_train: False,
                       attn_drop: 0.0,
                       ffd_drop: 0.0}

                fd = fd1
                fd.update(fd2)
                fd.update(fd3)
                loss_value_vl, acc_vl = sess.run([loss, accuracy],
                                                 feed_dict=fd)
                val_loss_avg += loss_value_vl
                val_acc_avg += acc_vl
                vl_step += 1
            # import pdb; pdb.set_trace()
            print('Epoch: {}, att_val: {}'.format(epoch, np.mean(att_val_train, axis=0)))
            print('Training: loss = %.5f, acc = %.5f | Val: loss = %.5f, acc = %.5f' %
                  (train_loss_avg / tr_step, train_acc_avg / tr_step,
                   val_loss_avg / vl_step, val_acc_avg / vl_step))

            if val_acc_avg / vl_step >= vacc_mx or val_loss_avg / vl_step <= vlss_mn:
                if val_acc_avg / vl_step >= vacc_mx and val_loss_avg / vl_step <= vlss_mn:
                    vacc_early_model = val_acc_avg / vl_step
                    vlss_early_model = val_loss_avg / vl_step
                    saver.save(sess, checkpt_file)
                vacc_mx = np.max((val_acc_avg / vl_step, vacc_mx))
                vlss_mn = np.min((val_loss_avg / vl_step, vlss_mn))
                curr_step = 0
            else:
                curr_step += 1
                if curr_step == patience:
                    print('Early stop! Min loss: ', vlss_mn,
                          ', Max accuracy: ', vacc_mx)
                    print('Early stop model validation loss: ',
                          vlss_early_model, ', accuracy: ', vacc_early_model)
                    break

            train_loss_avg = 0
            train_acc_avg = 0
            val_loss_avg = 0
            val_acc_avg = 0

        saver.restore(sess, checkpt_file)
        print('load model from : {}'.format(checkpt_file))
        ts_size = fea_list[0].shape[0]
        ts_step = 0
        ts_loss = 0.0
        ts_acc = 0.0

        while ts_step * batch_size < ts_size:
            # fd1 = {ftr_in: features[ts_step * batch_size:(ts_step + 1) * batch_size]}
            fd1 = {i: d[ts_step * batch_size:(ts_step + 1) * batch_size]
                   for i, d in zip(ftr_in_list, fea_list)}
            fd2 = {i: d[ts_step * batch_size:(ts_step + 1) * batch_size]
                   for i, d in zip(bias_in_list, biases_list)}
            fd3 = {lbl_in: y_test[ts_step * batch_size:(ts_step + 1) * batch_size],
                   msk_in: test_mask[ts_step * batch_size:(ts_step + 1) * batch_size],

                   is_train: False,
                   attn_drop: 0.0,
                   ffd_drop: 0.0}

            fd = fd1
            fd.update(fd2)
            fd.update(fd3)
            loss_value_ts, acc_ts, jhy_final_embedding = sess.run([loss, accuracy, final_embedding],
                                                                  feed_dict=fd)
            ts_loss += loss_value_ts
            ts_acc += acc_ts
            ts_step += 1

        print('****************************.....')
        print('Test loss:', ts_loss / ts_step,
              '; Test accuracy:', ts_acc / ts_step)

        #print('start knn, kmean.....')
        print('****************************.....')
        xx = np.expand_dims(jhy_final_embedding, axis=0)[test_mask]
        print(type(jhy_final_embedding))
        print(jhy_final_embedding.shape)
        print('****************************.....')


        # Specify the file path

        file_path_export = "/home/fatemeh/data/HCC/Export/Emb_g_5.txt" #12

        # Write the array to the file
        np.savetxt(file_path_export, jhy_final_embedding, fmt='%d')

        from numpy import linalg as LA

        # xx = xx / LA.norm(xx, axis=1)
        yy = y_test[test_mask]

        print('****************************.....')
        #print('xx: {}, yy: {}'.format(xx.shape, yy.shape))

        predictions = RClassification(jhy_final_embedding, labels)

        #predictions = SVMClassification(jhy_final_embedding, labels)

        #predictions = RFClassification(jhy_final_embedding, labels)

        #predictions = LRClassification(jhy_final_embedding, labels)

        #predictions = AdaBoostClassification(jhy_final_embedding, labels)

        #predictions = NaiveBayesClassification(jhy_final_embedding, labels)








        print('prediction:', predictions)
        #my_KNN(xx, yy)
        #my_Kmeans(xx, yy)

        sess.close()
