# -*- coding: UTF-8 -*-
"""
@Project ：计图 
@File ：main.py
@Author ：AnthonyZ
@Date ：2022/5/11 09:26
"""

import numpy as np
import jittor
import jittor.nn as nn

from jittor.dataset.mnist import MNIST
import jittor.transform as transform

import argparse
from tqdm import tqdm
import logging
import os

from PIL import Image


class Generator(nn.Module):
    def __init__(self, opt):
        super(Generator, self).__init__()
        self.args = opt
        self.embedding = nn.Embedding(10, 10)
        self.linear1 = nn.Linear((opt.latent_dim + opt.n_classes), 128)
        self.linear2 = nn.Linear(128, 256)
        self.batch_normal1 = nn.BatchNorm1d(256, 0.8)
        self.linear3 = nn.Linear(256, 512)
        self.batch_normal2 = nn.BatchNorm1d(512, 0.8)
        self.linear4 = nn.Linear(512, 1024)
        self.batch_normal3 = nn.BatchNorm1d(1024, 0.8)
        self.linear5 = nn.Linear(1024, int(opt.channels * opt.img_size * opt.img_size))
        self.relu = nn.LeakyReLU(0.2)
        self.tanh = nn.Tanh()

    def execute(self, noise, lab):
        x = jittor.contrib.concat((self.embedding(lab), noise), dim=1)
        x = self.linear1(x)
        x = self.linear2(self.relu(x))
        x = self.linear3(self.relu(self.batch_normal1(x)))
        x = self.linear4(self.relu(self.batch_normal2(x)))
        x = self.linear5(self.relu(self.batch_normal3(x)))
        x = self.tanh(x)
        x = x.view((-1, self.args.channels, self.args.img_size, self.args.img_size))
        return x


class Discriminator(nn.Module):
    def __init__(self, opt):
        super(Discriminator, self).__init__()
        self.args = opt
        self.embedding = nn.Embedding(10, 10)
        self.linear1 = nn.Linear((opt.n_classes + int(opt.channels * opt.img_size * opt.img_size)), 512)
        self.linear2 = nn.Linear(512, 512)
        self.linear3 = nn.Linear(512, 512)
        self.linear4 = nn.Linear(512, 1)
        self.drop_out = nn.Dropout(0.4)
        self.relu = nn.LeakyReLU(0.2)

    def execute(self, img, lab):
        img = img.view((img.shape[0], (- 1)))
        x = jittor.contrib.concat((img, self.embedding(lab)), dim=1)
        x = self.linear1(x)
        x = self.linear2(self.relu(x))
        x = self.linear3(self.drop_out(self.relu(x)))
        x = self.linear4(self.drop_out(self.relu(x)))
        return x

def init_logger():
    logger = logging.getLogger(__name__)
    logging.basicConfig(format='%(asctime)s - %(levelname)s - %(name)s -   %(message)s',
                        datefmt='%m/%d/%Y %H:%M:%S',
                        level=logging.INFO)

    logger.setLevel(logging.INFO)
    formatter = logging.Formatter(fmt='%(asctime)s - %(levelname)s - %(name)s -   %(message)s',
                                  datefmt='%m/%d/%Y %H:%M:%S')

    # 使用FileHandler输出到文件
    num_list = []
    for filename in os.listdir('./log'):
        if 'log' not in filename:
            continue
        num = int(filename.split('.')[0][3:])
        num_list.append(num)
    num = max(num_list) + 1 if num_list != [] else 1
    fh = logging.FileHandler('./log/log{}.txt'.format(num))
    fh.setLevel(logging.INFO)
    fh.setFormatter(formatter)

    logger.addHandler(fh)
    return logger

def save_image(img, path, nrow=10, padding=5):
    N, C, W, _ = img.shape
    if (N % nrow != 0):
        print("N %n row != 0")
        return
    ncol = int(N / nrow)
    img_all = []
    for i in range(ncol):
        img_ = []
        for j in range(nrow):
            img_.append(img[i * nrow + j])
            img_.append(np.zeros((C, W, padding)))
        img_all.append(np.concatenate(img_, 2))
        img_all.append(np.zeros((C, padding, img_all[0].shape[2])))
    img = np.concatenate(img_all, 1)
    img = np.concatenate([np.zeros((C, padding, img.shape[2])), img], 1)
    img = np.concatenate([np.zeros((C, img.shape[1] ,padding)), img], 2)
    min_ = img.min()
    max_ = img.max()
    img = (img - min_) / (max_ - min_) * 255
    img = img.transpose((1,2,0))
    if C == 3:
        img = img[: , : , : : -1]
    elif C == 1:
        img = img[: , : , 0]
    Image.fromarray(np.uint8(img)).save(path)

def train(args, train_loader):
    for i in range(args.epochs):
        train_tqdm = tqdm(enumerate(train_loader), desc="Epoch " + str(i))
        for _, (data, target) in train_tqdm:
            valid = jittor.ones([data.shape[0], 1]).float32().stop_grad()
            fake = jittor.zeros([data.shape[0]]).float32().stop_grad()

            real_imgs = jittor.array(data)
            labels = jittor.array(target)

            optimizer_G.zero_grad()

            noise = jittor.array(np.random.normal(0, 1, (data.shape[0], args.latent_dim))).float32()
            fake_labels = jittor.array(np.random.randint(0, args.n_classes, data.shape[0])).float32()

            fake_imgs = generator(noise, fake_labels)

            predict = discriminator(fake_imgs, fake_labels)
            loss1 = adversarial_loss(predict, valid)
            loss1.sync()
            optimizer_G.step(loss1)

            optimizer_D.zero_grad()

            predict_real = discriminator(real_imgs, labels)
            loss2 = adversarial_loss(predict_real, valid)
            predict_fake = discriminator(fake_imgs.stop_grad(), fake_labels)
            loss3 = adversarial_loss(predict_fake, fake)
            loss = (loss3 + loss2) / 2
            loss.sync()

            optimizer_D.step(loss)
        logger.info("[Epoch {}/{}] [D loss: {}] [G loss: {}]".format(i, args.epochs, loss.data, loss1.data))
        
        if i % 10 == 0:
            logger.info("Save example gen-image and model")
            generator.save("saved_models/generator_last.pkl")
            discriminator.save("saved_models/discriminator_last.pkl")
            n_row = 10
            batches_done=i
            labels_temp = jittor.array(np.array([num for _ in range(n_row) for num in range(n_row)])).float32().stop_grad()
            gen_imgs = generator(jittor.array(np.random.normal(0, 1, (n_row ** 2, opt.latent_dim))).float32().stop_grad(), labels_temp)
            save_image(gen_imgs.numpy(), "./example/%d.png" % batches_done, nrow=n_row)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--epochs', type=int, default=100, help='number of epochs of training')
    parser.add_argument('--batch_size', type=int, default=64, help='size of the batches')
    parser.add_argument('--lr', type=float, default=0.0002, help='adam: learning rate')
    parser.add_argument('--latent_dim', type=int, default=100, help='dimensionality of the latent space')
    parser.add_argument('--n_classes', type=int, default=10, help='number of classes for dataset')
    parser.add_argument('--img_size', type=int, default=32, help='size of each image dimension')
    parser.add_argument('--channels', type=int, default=1, help='number of image channels')
    parser.add_argument('--device', type=str, default='CPU', help='Training device')
    opt = parser.parse_args()

    if opt.device == 'GPU':
        jittor.flags.use_cuda = 1

    os.makedirs("example", exist_ok=True)
    os.makedirs("log", exist_ok=True)
    os.makedirs("saved_models", exist_ok=True)

    number = '18713012939'
    z = jittor.array(np.random.normal(0, 1, (len(number), opt.latent_dim))).float32().stop_grad()
    labels = jittor.array(np.array([int(number[num]) for num in range(len(number))])).float32().stop_grad()

    transform = transform.Compose([
        transform.Resize(opt.img_size),
        transform.Gray(),
        transform.ImageNormalize(mean=[0.5], std=[0.5]),
    ])
    dataloader = MNIST(train=True, transform=transform).set_attrs(batch_size=opt.batch_size, shuffle=True)
    
    generator = Generator(opt)
    discriminator = Discriminator(opt)

    optimizer_G = nn.Adam(generator.parameters(), lr=opt.lr, betas=(0.5, 0.999))
    optimizer_D = nn.Adam(discriminator.parameters(), lr=opt.lr, betas=(0.5, 0.999))

    adversarial_loss = jittor.nn.MSELoss()

    logger = init_logger()

    train(opt, dataloader)

    generator.load('saved_models/generator_last.pkl')
    discriminator.load('saved_models/discriminator_last.pkl')

    gen_imgs = generator(jittor.array(np.random.normal(0, 1, (len(number), opt.latent_dim))).float32().stop_grad(), jittor.array(np.array([int(number[num]) for num in range(len(number))])).float32().stop_grad())
    img_array = gen_imgs.data.transpose((1, 2, 0, 3))[0].reshape((gen_imgs.shape[2], -1))
    img_array = (img_array - img_array.min()) / (img_array.max() - img_array.min()) * 255
    Image.fromarray(np.uint8(img_array)).save("result.png")
    

