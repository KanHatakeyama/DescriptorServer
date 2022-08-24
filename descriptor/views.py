from django.shortcuts import render
from django.http import HttpResponse
# Create your views here.


def calc(request):
    smiles = (request.GET["SMILES"])
    ret = smiles
    return HttpResponse(ret)
