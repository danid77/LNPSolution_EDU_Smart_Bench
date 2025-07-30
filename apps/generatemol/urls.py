from django.urls import path

from . import views

app_name = 'generatemol'

urlpatterns = [
    path('molsparkInputPage/', views.molsparkInputPage, name='molsparkInputPage'),
    path('molsparkOutputPage/', views.molsparkOutputPage, name='molsparkOutputPage'),
    path('molsparkProcess/', views.molsparkProcess, name='molsparkProcess'),
    
    path('chemDrawInputPage/', views.chemDrawInputPage, name='chemDrawInputPage'),
    path('chemDrawProcess/', views.chemDrawProcess, name='chemDrawProcess'),
    path('chemDrawOutputPage/', views.chemDrawOutputPage, name='chemDrawOutputPage'),
    path('chemDrawSmilesTo3d/', views.chemDrawSmilesTo3d, name='chemDrawSmilesTo3d'),
    
    path('pocket2MolInputPage/', views.pocket2MolInputPage, name='pocket2MolInputPage'),
    path('pocket2MolProcess/', views.pocket2MolProcess, name='pocket2MolProcess'),
    path('pocket2MolOutputPage/', views.pocket2MolOutputPage, name='pocket2MolOutputPage'),
]