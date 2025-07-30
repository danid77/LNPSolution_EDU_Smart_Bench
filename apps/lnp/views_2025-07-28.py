from django.shortcuts import render
from django.urls import reverse
from django.http import JsonResponse
from django.http import HttpResponse
import json, os

# Create your views here.

def home(request):
    return render(request, 'lnp/home.html')

def singleModule(request):
    cards = [
        {
            "title": "Protein Structure Prediction - Monomer",
            "icon": "img/single/protein_gen2.png",
            "color": "#198754",
            "url": reverse('lnp:proteinMonomerHome'),
            "delay": "0.1s"
        },
        {
            "title": "Chem draw : ketcher",
            "icon": "img/nvidia_molmim.jpg",
            "color": "#198754",
            "url": reverse('generatemol:chemDrawInputPage'),
            "delay": "0.2s"
        },
        {
            "title": "Complex Prediction / Docking",
            "icon": "img/single/complex2.png",
            "color": "#198754",
            "url": reverse('lnp:complexHome'),
            "delay": "0.3s"
        },
        {
            "title": "Protein Structure Prediction - Multimer",
            "icon": "img/single/protein_gen2.png",
            "color": "#198754",
            "url": reverse('lnp:proteinMultimerHome'),
            "delay": "0.4s"
        },
        {
            "title": "Small Molecule Generation / Optimization",
            "icon": "img/single/molecule.png",
            "color": "#198754",
            "url": reverse('lnp:moleculeHome'),
            "delay": "0.5s"
        },
        {
            "title": "Protein generation - ProteinMPNN",
            "icon": "img/nvidia_proteinmpnn.webp",
            "color": "#198754",
            "url": reverse('nim:nimProteinmpnnInputPage'),
            "delay": "0.6s"
        },
        {
            "title": "Pocket-based compound generation",
            "icon": "img/nvidia_molmim.jpg",
            "color": "#198754",
            "url": reverse('lnp:test'),
            "delay": "0.7s"
        },
        {
            "title": "From Binding to Protein - Rfdiffusion",
            "icon": "img/nvidia_rfdiffusion.webp",
            "color": "#198754",
            "url": reverse('nim:nimRfdiffusionInputPage'),
            "delay": "0.8s"
        },
        # {
        #     "title": "Peptide discovery",
        #     "icon": "img/single/peptide2.png",
        #     "color": "#198754",
        #     "url": reverse('lnp:peptideHome'),
        #     "delay": "0.1s"
        # },
        # {
        #     "title": "Interaction Analysis & Visualization",
        #     "icon": "img/single/interaction.png",
        #     "color": "#198754",
        #     "url": reverse('lnp:interactionHome'),
        #     "delay": "0.1s"
        # },
        # {
        #     "title": "Dimensionality Reduction & Visualization",
        #     "icon": "img/single/umap_img.png",
        #     "color": "#198754",
        #     "url": reverse('lnp:dimensionalityHome'),
        #     "delay": "0.1s"
        # },
    ]
    return render(request, 'lnp/single_module.html', {'cards': cards})

def workFlow(request):
    return render(request, 'lnp/work_flow.html')

def proteinMonomerHome(request):
    cards = [
        {
            "title": "Alphafold2",
            "icon": "img/nvidia_alphafold2.jpg",
            "color": "#198754",
            "url": reverse('nim:nimAlphafoldInputPage'),
            "delay": "0.1s"
        },
        {
            "title": "Openfold",
            "icon": "img/nvidia_openfold2.webp",
            "color": "#198754",
            "url": reverse('nim:nimOpenfoldInputPageSample'),
            "delay": "0.5s"
        },
        {
            "title": "ESMfold",
            "icon": "img/nvidia_esmfold.webp",
            "color": "#198754",
            "url": reverse('nim:nimEsmfoldInputPage'),
            "delay": "0.7s"
        },
        # {
        #     "title": "Boltz",
        #     "icon": "img/single_module/boltz.png",
        #     "color": "#4b9153",
        #     "url": reverse('protein:boltzInputPage'),
        #     "delay": "0.9s"
        # },
        {
            "title": "Boltz2",
            "icon": "img/single_module/boltz.png",
            "color": "#4b9153",
            "url": reverse('protein:boltz2InputPage'),
            "delay": "0.9s"
        },
        # {
        #     "title": "Chai - 1",
        #     "icon": "img/single_module/chai.png",
        #     "color": "#4b9153",
        #     "url": reverse('protein:chaiInputPage'),
        #     "delay": "1.1s"
        # },
    ]
    return render(request, 'lnp/category.html', {'cards': cards})

def proteinMultimerHome(request):
    cards = [
        {
            "title": "Alphafold2-Multimer",
            "icon": "img/nvidia_alphafold2_multimer.jpg",
            "color": "#198754",
            "url": reverse('nim:nimAlphafoldMultiInputPage'),
            "delay": "0.3s"
        },
        {
            "title": "Boltz2-Multimer",
            "icon": "img/single_module/boltz.png",
            "color": "#4b9153",
            "url": reverse('protein:boltz2MultiInputPage'),
            "delay": "0.9s"
        },
    ]
    return render(request, 'lnp/category.html', {'cards': cards})

def moleculeHome(request):
    cards = [
        {
            "title": "Chem draw : ketcher",
            "icon": "img/nvidia_molmim.jpg",
            "color": "#198754",
            "url": reverse('generatemol:chemDrawInputPage'),
            "delay": "0.1s"
        },
        {
            "title": "Molmin",
            "icon": "img/nvidia_molmim.jpg",
            "color": "#198754",
            "url": reverse('nim:nimMolminInputPageSample'),
            "delay": "0.1s"
        },
        {
            "title": "Genmol",
            "icon": "img/nvidia_genmol.webp",
            "color": "#198754",
            "url": reverse('nim:nimGenmolInputPageSample'),
            "delay": "0.3s"
        },
                {
            "title": "MolSpark",
            "icon": "img/single/molecule.png",
            "color": "#dc3545",
            "url": reverse('generatemol:molsparkInputPage'),
            "delay": "0.5s"
        },
        # {
        #     "title": "MolSpark",
        #     "icon": "img/single/molecule.png",
        #     "color": "#dc3545",
        #     "url": "http://lnpsolution.iptime.org:8601/",
        #     "delay": "0.5s"
        # },
    ]
    return render(request, 'lnp/category.html', {'cards': cards})

def complexHome(request):
    cards = [
        {
            "title": "Diffdock",
            "icon": "img/nvidia_diffdock.jpg",
            "color": "#198754",
            "url": reverse('nim:nimDiffdockInputPage'),
            "delay": "0.1s"
        },
        # {
        #     "title": "Boltz",
        #     "icon": "img/single_module/boltz_complex.png",
        #     "color": "#4b9153",
        #     "url": reverse('protein:boltzComplexInputPage'),
        #     "delay": "0.3s"
        # },
        {
            "title": "Boltz2",
            "icon": "img/single_module/boltz_complex.png",
            "color": "#4b9153",
            "url": reverse('protein:boltz2ComplexInputPage'),
            "delay": "0.3s"
        },
        # {
        #     "title": "Chai - 1",
        #     "icon": "img/single_module/chai_complex.png",
        #     "color": "#4b9153",
        #     "url": reverse('protein:chaiComplexInputPage'),
        #     "delay": "0.5s"
        # },
        # {
        #     "title": "GNINA",
        #     "icon": "img/nvidia_molmim.jpg",
        #     "color": "#198754",
        #     "url": reverse('lnp:test'),
        #     "delay": "0.7s"
        # },
    ]
    return render(request, 'lnp/category.html', {'cards': cards})

def peptideHome(request):
    cards = [
        {
            "title": "Peptide generation - ProteinMPNN",
            "icon": "img/nvidia_proteinmpnn.webp",
            "color": "#198754",
            "url": reverse('nim:nimProteinmpnnInputPage'),
            "delay": "0.1s"
        },
        {
            "title": "Rfdiffusion",
            "icon": "img/nvidia_rfdiffusion.webp",
            "color": "#198754",
            "url": reverse('nim:nimRfdiffusionInputPage'),
            "delay": "0.3s"
        },
        # {
        #     "title": "HydraAMP",
        #     "icon": "img/single/peptide2.png",
        #     "color": "#4b9153",
        #     "url": "http://lnpsolution.iptime.org:8610/",
        #     "delay": "0.5s"
        # },
    ]
    return render(request, 'lnp/category.html', {'cards': cards})

# def interactionHome(request):
#     cards = [
#         {
#             "title": "PLIP",
#             "icon": "img/logo3.png",
#             "color": "#0d6efd",
#             "url": reverse('plip:plipInputPage'),
#             "delay": "0.1s"
#         },
#         {
#             "title": "Visualizing: Protein & Ligand Interaction",
#             "icon": "img/logo3.png",
#             "color": "#0d6efd",
#             "url": reverse('visualzing:visualzingInputPage'),
#             "delay": "0.2s"
#         },
#     ]
#     return render(request, 'lnp/category.html', {'cards': cards})

def interactionHome(request):
    cards = [
        {
            "title": "Visualizing",
            "icon": "img/single/interaction_3.png",
            "color": "#0d6efd",
            "url": reverse('visualzing:visualzingInputPage'),
            "delay": "0.1s"
        },
        {
            "title": "PLIP",
            "icon": "img/single_module/plip_3.png",
            "color": "#0d6efd",
            "url": reverse('plip:plipInputPage'),
            "delay": "0.3s"
        },
    ]
    return render(request, 'lnp/category.html', {'cards': cards})

def dimensionalityHome(request):
    cards = [
        {
            "title": "Molecular Data Visualization",
            "icon": "img/single/umap_img.png",
            "color": "#4b9153",
            "url": "http://lnpsolution.iptime.org:8502/",
            "delay": "0.1s"
        },
    ]
    return render(request, 'lnp/category.html', {'cards': cards})


def test(request):
    return render(request, 'lnp/test.html')


# def lnpWork(request):
#     cards = [
#         {
#             "title": "Chemical Converter",
#             "icon": "bi bi-capsule",
#             "color": "#4b9153",
#             "link": "http://lnplicense.iptime.org:8606/",
#             "delay": "0.5s"
#         },
#     ]

#     return render(request, 'lnp/lnp_work.html', {'cards': cards})

# def mix(request):
#     cards = [
#         {
#             "title": "protein",
#             "icon": "bi bi-cpu",
#             "color": "#198754",
#             "url": reverse('protein:proteinInputPage'),
#             "delay": "0.1s"
#         },
#         {
#             "title": "Complex Analysis",
#             "icon": "bi bi-bar-chart-line",
#             "color": "#0d6efd",
#             "url": reverse('screening:complexInputPage'),
#             "delay": "0.3s"
#         },
#         {
#             "title": "Chemical generator & analysis",
#             "icon": "bi bi-stars",
#             "color": "#0d6efd",
#             "url": reverse('umap:umapInputPageSample'),
#             "delay": "0.1s"
#         },
        
#     ]

#     return render(request, 'lnp/mix.html', {'cards': cards})