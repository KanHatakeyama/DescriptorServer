from django.shortcuts import render
from django.http import HttpResponse
from django.http import JsonResponse
from .calculators.wrapper import fetch_descriptor, process_smiles
from django.views import View
import pandas as pd
import json
# Create your views here.


class PostView(View):
    def get(self, request, *args, **kwargs):
        context = {
            'SMILES': "C\nCOC",
        }
        return render(request, 'post.html', context)

    def post(self, request, *args, **kwargs):
        smiles_list = request.POST['SMILES']

        if "opt[]" in request.POST:
            option_list = request.POST.getlist("opt[]")
        else:
            option_list = []

        processed_dict = process_smiles(smiles_list, option_list)

        # export json
        if "json" in request.POST:
            return JsonResponse((processed_dict))

        # export csv
        df = pd.DataFrame.from_dict(processed_dict).T
        filename = request.POST["filename"]

        response = HttpResponse(content_type='text/csv; charset=utf8')
        response['Content-Disposition'] = f'attachment; filename={filename}'
        df.to_csv(path_or_buf=response, encoding='utf_8_sig', index=None)
        return response


post_view = PostView.as_view()
