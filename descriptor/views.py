from django.shortcuts import render
from django.http import HttpResponse
from django.http import JsonResponse
from .calculators.wrapper import fetch_descriptor
# Create your views here.


def calc(request):
    ret = fetch_descriptor(request)

    return JsonResponse(ret)
