import socket
import sys

def eval_srepr(srepr):
	return eval(srepr);

def gethostbyname(hostname):
	return socket.gethostbyname(hostname);
