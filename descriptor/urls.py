from django.urls import path

from . import views

urlpatterns = [
    path('calc', views.calc),
    path('post', views.post_view),
]
