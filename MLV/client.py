# to send request to MLV server

import json
from socket import *
from select import select


HOST = '127.0.0.1'
PORT = 35520
BUFSIZE = 1024
ADDR = (HOST, PORT)

clientSocket = socket(AF_INET, SOCK_STREAM)
try:
  clientSocket.connect(ADDR)
except Exception as e:
  raise Exception('Cannot connect to MLV task server!!')




def send(name):
  msg = json.dumps({'name': name})
  clientSocket.send(msg)

def close():
  clientSocket.close()
