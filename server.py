# real processes doing RWR task

from socket import *
from select import *
import threading, time
from MLV import MLV
from time import ctime
import json

tasks = []
bRunning = True

def task_scanner():
  print '# start: task scanner'
  while bRunning:
    if len(tasks):
      task = tasks[0]
      tasks.pop(0)
      MLV.generate_graph(task['name'])
      s = MLV.load_state(task['name'])
      s['result'] = 1
      MLV.save_state(s)

def start_socket_server():
  HOST = ''
  PORT = 35520
  BUFSIZE = 1024
  ADDR = (HOST, PORT)
  serverSocket = socket(AF_INET, SOCK_STREAM)
  serverSocket.bind(ADDR)
  serverSocket.listen(5)
  connection_list = [serverSocket]
  global tasks
  while True:
    try:
      print '# start: socket server'
      read_socket, write_socket, error_socket = select(connection_list, [], [], 5)
      for sock in read_socket:
        if sock == serverSocket:
          clientSocket, addr_info = serverSocket.accept()
          connection_list.append(clientSocket)
        else:
          data = sock.recv(BUFSIZE)
          if (not data):
            connection_list.remove(sock)
            sock.close()
            continue
          task = json.loads(data)
          tasks.append(task)
    except KeyboardInterrupt:
      serverSocket.close()

if (__name__=="__main__"):
  # start task scanner first
  t1 = threading.Thread(target=task_scanner)
  t1.daemon = True
  t1.start()
  # open socket server (till end)
  start_socket_server()
  # end of server; exit thread
  bRunning = False
  print 'bye'
