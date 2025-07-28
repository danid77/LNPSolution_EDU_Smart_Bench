from django.urls import path

from . import views

app_name = 'protein'

urlpatterns = [
    path('proteinInputPage/', views.proteinInputPageSample, name='proteinInputPageSample'),
    path('proteinProcess/', views.proteinProcessSample, name='proteinProcessSample'),
    path('proteinOutputPage/', views.proteinOutputPage, name='proteinOutputPage'),
    
    path('boltzInputPage/', views.boltzInputPage, name='boltzInputPage'),
    path('boltzProcess/', views.boltzProcess, name='boltzProcess'),
    path('boltzOutputPage/', views.boltzOutputPage, name='boltzOutputPage'),
    
    path('boltz2OutputListPage/', views.boltz2OutputListPage, name='boltz2OutputListPage'),
    path('boltz2InputPage/', views.boltz2InputPage, name='boltz2InputPage'),
    path('boltz2Process/', views.boltz2Process, name='boltz2Process'),
    path('boltz2OutputPage/', views.boltz2OutputPage, name='boltz2OutputPage'),
    
    path('boltz2MultiInputPage/', views.boltz2MultiInputPage, name='boltz2MultiInputPage'),
    path('boltz2MultiProcess/', views.boltz2MultiProcess, name='boltz2MultiProcess'),
    path('boltz2MultiOutputPage/', views.boltz2MultiOutputPage, name='boltz2MultiOutputPage'),
    
    path('chaiInputPage/', views.chaiInputPage, name='chaiInputPage'),
    path('chaiProcess/', views.chaiProcess, name='chaiProcess'),
    path('chaiOutputPage/', views.chaiOutputPage, name='chaiOutputPage'),
    
    path('boltzComplexInputPage/', views.boltzComplexInputPage, name='boltzComplexInputPage'),
    path('boltzComplexProcess/', views.boltzComplexProcess, name='boltzComplexProcess'),
    path('boltzComplexOutputPage/', views.boltzComplexOutputPage, name='boltzComplexOutputPage'),
    
    path('boltz2ComplexInputPage/', views.boltz2ComplexInputPage, name='boltz2ComplexInputPage'),
    path('boltz2ComplexProcess/', views.boltz2ComplexProcess, name='boltz2ComplexProcess'),
    path('boltz2ComplexOutputPage/', views.boltz2ComplexOutputPage, name='boltz2ComplexOutputPage'),
    
    path('chaiComplexInputPage/', views.chaiComplexInputPage, name='chaiComplexInputPage'),
    path('chaiComplexProcess/', views.chaiComplexProcess, name='chaiComplexProcess'),
    path('chaiComplexOutputPage/', views.chaiComplexOutputPage, name='chaiComplexOutputPage'),
]