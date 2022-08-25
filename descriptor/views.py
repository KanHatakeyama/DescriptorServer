from django.shortcuts import render
from django.http import HttpResponse
from django.http import JsonResponse
from .calculators.wrapper import fetch_descriptor
from django.views import View
# Create your views here.


def calc(request):
    ret = fetch_descriptor(request)

    return JsonResponse(ret)


class PostView(View):
    def get(self, request, *args, **kwargs):
        context = {
            'SMILES': "C\nCOC",
        }
        return render(request, 'post.html', context)

    def post(self, request, *args, **kwargs):
        smiles_list = request.POST['SMILES']
        option_list = request.POST["opt"]
        print(smiles_list, option_list)

        context = {
            'SMILES': request.POST['SMILES']+"a",
        }

        return render(request, 'post.html', context)


post_view = PostView.as_view()
