from django.shortcuts import render
from django.http import JsonResponse
from django.urls import reverse
from django.utils.http import urlsafe_base64_decode
import requests
from pathlib import Path
import json, os
import pandas as pd
import numpy as np
import tempfile
from datetime import datetime
import subprocess
import time
from urllib.parse import unquote

from apps import getKey, outputPath, outputSamplePath, sampleList, multiSampleList
from apps.umap.services import molMin, genMol, molsparkProcess
# Create your views here.
        
def nimOpenfoldInputPageSample(request):
    databases = ['Uniref30_2302', 'PDB70_220313', 'colabfold_envdb_202108']
    return render(request, 'nim/openfold/input_sample.html', {'databases': databases, 'samples_json': json.dumps(sampleList())})

def nimOpenfoldProcessSample(request):
    if request.method == "POST":    
        try:
            sample_dir = f"{outputSamplePath()}/nim/openfold_sample"
            
            data = json.loads(request.body)
            sample_name = data.get('sampleName', '')  
            # print(sample_name)
            time.sleep(1)  # 5초 동안 멈춤
            results_dir = f"{sample_dir}/{sample_name}"
            
            return JsonResponse({"status": "done", "results_dir" : results_dir})  # ← 클라이언트에 저장 완료만 알림
        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)
        
def nimOpenfoldList(request):   
    try:
        results_dir = request.GET.get('results_dir')
        
        with open(f"{results_dir}/openfold2.json", 'r') as f:
            openfold = json.load(f)
        if openfold is None:
            return JsonResponse({"error": "세션에 저장된 결과가 없습니다."}, status=404)
        
        cards = []
        for i in range(len(openfold['structures_in_ranked_order'])):
            title_pretty = f"Model {i + 1}"
            url = f"{reverse('nim:nimOpenfoldOutputPage')}?results_dir={results_dir}&num={i}"
                
            cards.append({
                "title": title_pretty,
                "icon": "img/nvidia_openfold2.webp",
                "color": "#198754",
                "url": url,
                "delay": f"{0.1 + i * 0.1:.1f}s"
            })
        
        return render(request, 'list.html', {"cards": cards, "tool" : "Openfold2"})
    except Exception as e:
        print("에러 발생:", e)
        return JsonResponse({"error": str(e)}, status=500)
    
# def boltz2OutputListPage(request):
#     results_dir = request.GET.get('results_dir')
#     source = request.GET.get('source')

#     pdb_paths = glob.glob(f"{results_dir}/*/predictions/*/*.pdb")
#     pdb_paths.sort()

#     if not pdb_paths:
#         return render(request, 'boltz2/list.html', {
#             "cards": [],
#             "pdb_paths": [],
#             "error_message": "❌ 생성된 모델이 없습니다. 입력값을 확인하거나 다시 시도해 주세요."
#         })

#     cards = []
#     for i, pdb_path in enumerate(pdb_paths):
#         file_name = os.path.basename(pdb_path)
#         title = file_name.replace(".pdb", "")
#         title_pretty = f"Model {i + 1}"
#         rel_pdb_path = pdb_path.replace("/opt/platform/LNPSolution_AI_Smart_Bench", "")

#         if source == "single" :
#             url = f"{reverse('protein:boltz2OutputPage')}?pdb_path={rel_pdb_path}"
#         elif source == "multi":
#             url = f"{reverse('protein:boltz2MultiOutputPage')}?pdb_path={rel_pdb_path}"
#         elif source == "complex_single":
#             url = f"{reverse('protein:boltz2ComplexOutputPage')}?pdb_path={rel_pdb_path}"
#         elif source == "complex_multi":
#             url = f"{reverse('protein:boltz2MultiComplexOutputPage')}?pdb_path={rel_pdb_path}"
            
#         cards.append({
#             "title": title_pretty,
#             "icon": "img/single_module/boltz.png",
#             "color": "#198754",
#             "url": url,
#             "delay": f"{0.1 + i * 0.1:.1f}s"
#         })

#     return render(request, 'boltz2/list.html', {
#         "cards": cards,
#         "pdb_paths": pdb_paths
#     })

def nimOpenfoldOutputPage(request):
    results_dir = request.GET.get('results_dir')
    num = int(request.GET.get('num'))
    
    with open(f"{results_dir}/openfold2.json", 'r') as f:
        openfold = json.load(f)
    if openfold is None:
        return JsonResponse({"error": "세션에 저장된 결과가 없습니다."}, status=404)
    
    result_pdb = openfold['structures_in_ranked_order'][num]['structure']
    confidence = openfold['structures_in_ranked_order'][num]['confidence']
    return render(request, 'nim/openfold/result.html', {'result_pdb' : result_pdb, "confidence" : confidence})


def nimMsaInputPage(request):
    sections = ["OpenFold2", "MMseqs2 MSA", "MolMin", "DiffDock Inputs"]
    databases = ["Uniref30_2303", "PDB70_230113", "colabfold_envdb_202108"]
    return render(request, 'nim/msa/input.html', {
        'sections': sections,
        'databases': databases
    })


def nimAlphafoldInputPage(request):
    databases = ["uniref90", "small_bfd", "mgnify"]
    return render(request, 'nim/af2/input.html', {'databases': databases, 'samples_json': json.dumps(sampleList())})
    
def nimAlphafoldProcess(request):
    if request.method == "POST":    
        try:
            sample_dir = f"{outputSamplePath()}/nim/af2"
            
            data = json.loads(request.body)
            sample_name = data.get('sampleName', '')  
            # print(sample_name)
            time.sleep(1)  # 5초 동안 멈춤
            results_dir = f"{sample_dir}/{sample_name}"
            
            return JsonResponse({"status": "done", "results_dir" : results_dir})  # ← 클라이언트에 저장 완료만 알림
        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)
        
def nimAlphafoldList(request):   
    try:
        results_dir = request.GET.get('results_dir')
        
        with open(f"{results_dir}/alphafold2.json", 'r') as f:
            fold = json.load(f)
        if fold is None:
            return JsonResponse({"error": "세션에 저장된 결과가 없습니다."}, status=404)
        
        cards = []
        for i in range(len(fold)):
            title_pretty = f"Model {i + 1}"
            url = f"{reverse('nim:nimAlphafoldOutputPage')}?results_dir={results_dir}&num={i}"
                
            cards.append({
                "title": title_pretty,
                "icon": "img/nvidia_alphafold2.jpg",
                "color": "#198754",
                "url": url,
                "delay": f"{0.1 + i * 0.1:.1f}s"
            })
        
        return render(request, 'list.html', {"cards": cards, "tool" : "Alphafold2"})
    except Exception as e:
        print("에러 발생:", e)
        return JsonResponse({"error": str(e)}, status=500)        

def nimAlphafoldOutputPage(request):
    results_dir = request.GET.get('results_dir')
    num = int(request.GET.get('num'))
    
    with open(f"{results_dir}/alphafold2.json", 'r') as f:
        fold = json.load(f)
    if fold is None:
        return JsonResponse({"error": "세션에 저장된 결과가 없습니다."}, status=404)

    return render(request, 'nim/af2/result.html', {"result_pdb" : fold[num], "num" : num + 1})


# ##############################################################################################
def nimEsmfoldInputPage(request):
    return render(request, 'nim/esmfold/input.html', {'samples_json': json.dumps(sampleList())})


def nimEsmfoldProcess(request):
    if request.method == "POST":    
        try:
            sample_dir = f"{outputSamplePath()}/nim/esmfold"
            
            data = json.loads(request.body)
            sample_name = data.get('sampleName', '')  
            # print(sample_name)
            time.sleep(1)  # 5초 동안 멈춤
            results_dir = f"{sample_dir}/{sample_name}"
            
            return JsonResponse({"status": "done", "results_dir" : results_dir})  # ← 클라이언트에 저장 완료만 알림
        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)

def nimEsmfoldList(request):   
    try:
        results_dir = request.GET.get('results_dir')
        
        with open(f"{results_dir}/esmfold.json", 'r') as f:
            fold = json.load(f)
        if fold is None:
            return JsonResponse({"error": "세션에 저장된 결과가 없습니다."}, status=404)
        
        cards = []
        for i in range(len(fold['structures'])):
            title_pretty = f"Model {i + 1}"
            url = f"{reverse('nim:nimEsmOutputPage')}?results_dir={results_dir}&num={i}"
                
            cards.append({
                "title": title_pretty,
                "icon": "img/nvidia_alphafold2.jpg",
                "color": "#198754",
                "url": url,
                "delay": f"{0.1 + i * 0.1:.1f}s"
            })
        
        return render(request, 'list.html', {"cards": cards, "tool" : "ESMfold"})
    except Exception as e:
        print("에러 발생:", e)
        return JsonResponse({"error": str(e)}, status=500)  

def nimEsmOutputPage(request):
    results_dir = request.GET.get('results_dir')
    num = int(request.GET.get('num'))
    
    with open(f"{results_dir}/esmfold.json", 'r') as f:
        fold = json.load(f)
    if fold is None:
        return JsonResponse({"error": "세션에 저장된 결과가 없습니다."}, status=404)

    return render(request, 'nim/af2/result.html', {"result_pdb" : fold['structures'][num]['pdb'], "num" : num + 1})



def nimGenmolInputPageSample(request):
    return render(request, 'nim/genmol/input_sample.html')

def nimGenmolProcessSample(request):
    if request.method == "POST":    
        try:
            sample_dir = f"{outputPath()}/nim/genmol_sample"
            data = json.loads(request.body)
            sample_name = data.get('molSeq', '')
            time.sleep(1)  # 5초 동안 멈춤
            results_dir = f"{sample_dir}/{sample_name}"
            return JsonResponse({"status": "done", "results_dir" : results_dir})  # ← 클라이언트에 저장 완료만 알림
        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)

def nimGenmolOutputPage(request):
    results_dir = request.GET.get('results_dir')
    if not results_dir or not os.path.isdir(results_dir):
        return JsonResponse({'error': 'Invalid results_dir'}, status=400)
    
    genmol_df = pd.read_csv(f"{results_dir}/genmol_result.csv")
    genmol_cal_df = pd.read_csv(f"{results_dir}/Rdkit_calculate_result.csv").drop(columns=[f"fp_{i+1}" for i in range(2048)])
    return render(request, 'nim/genmol/result.html', {
                                                      "genmol_num" : len(genmol_df), 
                                                      "smiles" : genmol_df['smiles'].tolist(), 
                                                      "cal_data": genmol_cal_df.replace({np.nan: None}).to_dict(orient="records")
                                                      })



        
def nimMolminInputPageSample(request):
    return render(request, 'nim/molmin/input_sample.html')


def nimMolminProcessSample(request):
    if request.method == "POST":    
        try:
            sample_dir = f"{outputSamplePath()}/nim/molmin_sample"
            data = json.loads(request.body)
            sample_name = data.get('molSeq', '')
            time.sleep(1)  # 5초 동안 멈춤
            results_dir = f"{sample_dir}/{sample_name}"
            return JsonResponse({"status": "done", "results_dir" : results_dir})  # ← 클라이언트에 저장 완료만 알림
        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)

def nimMolminOutputPage(request):
    results_dir = request.GET.get('results_dir')
    if not results_dir or not os.path.isdir(results_dir):
        return JsonResponse({'error': 'Invalid results_dir'}, status=400)
    
    result_df = pd.read_csv(f"{results_dir}/molmin_result.csv")
    result_cal_df = pd.read_csv(f"{results_dir}/Rdkit_calculate_result.csv").drop(columns=[f"fp_{i+1}" for i in range(2048)])
    return render(request, 'nim/molmin/result.html', {"num" : len(result_df), 
                                                      "smiles" : result_df['smiles'].tolist(),
                                                      "cal_data": result_cal_df.replace({np.nan: None}).to_dict(orient="records")                              
                                                      })

# def nimAlphafoldMultiInputPage(request):
#     databases = ["uniref90", "small_bfd", "mgnify"]
    
#     sequences1 = ["MGSKKLKRVGLSQELCDRLSRHQILTCQDFLCLSPLELMKVTGLSYRGVHELLCMVSRACAPKMQTAYGIKAQRSADFSPAFLSTTLSALDEALHGGVACGSLTEITGPPGCGKTQFCIMMSILATLPTNMGGLEGAVVYIDTESAFSAERLVEIAESRFPRYFNTEEKLLLTSSKVHLYRELTCDEVLQRIESLEEEIISKGIKLVILDSVASVVRKEFDAQLQGNLKERNKFLAREASSLKYLAEEFSIPVILTNQITTHLSGALASQADLVSPADDLSLSEGTSGSSCVIAALGNTWSHSVNTRLILQYLDSERRQILIAKSPLAPFTSFVYTIKEEGLVLQAYGNS", 
#               "GSHMAQPRPPFHITIPIYPGVDLLDVAAPVELFSWMADAWKARATTITLAAEHLTPLKTRDGLTLTPQRQFADYADAAAPQPQTHLLWVPGGAPDVLRKLMRGGPYLDFLKAQSAGADHVSSVCEGALLLAAAGLLDGYRATTHWAFIPCLQQFPAIKVAEGFPRYVIDGNRITGGGISSGLAEALAIVARVAGQDIAKHVQMITQYFPDPPFEQTIVPATHCPLQA"]  
#     sequences2 = ["QSALTQPASVSGSPGQSITISCTGTSSDVGGYNYVSWYQQHPGKAPKLMIYEVTNRPSGVSNRFSGSRSGNTASLTISGLQAEDEADYYCSSYTSSSLYVFGTGTKVAVLGQPKANPTVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADSSPVKAGVETTKPSKQSNNKYAASSYLSLTPEQWKSHRSYSCQVTHEGSTVEKTVAPTECS", 
#                 "QVHLVQSGAEVKKPGSSVKVSCKASGGTFSSCAISWVRQAPGQGLEWMGRIIPILGIANYAQKFQGRVTITADKSTSTAYMELSSLRSEDTAVYYCARGWEFGSGSYYRTDYYYYAMDVWGQGTTVTVSSASTKGPSVFPLAPCSRSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPK",
#                 "SNMNNTASWFTALTQHGKEDLKFPRGQGVPINTNSSPDDQIGYYRRATRRIRGGDGKMKDLSPRWYFYYLGTGPEAGLPYGANKDGIIWVATEGALNTPKDHIGTRNPANNAAIVLQLPQGTTLPKGFYA"]  
#     sequences3 = ["MGSKKLKRVGLSQELCDRLSRHQILTCQDFLCLSPLELMKVTGLSYRGVHELLCMVSRACAPKMQTAYGIKAQRSADFSPAFLSTTLSALDEALHGGVACGSLTEITGPPGCGKTQFCIMMSILATLPTNMGGLEGAVVYIDTESAFSAERLVEIAESRFPRYFNTEEKLLLTSSKVHLYRELTCDEVLQRIESLEEEIISKGIKLVILDSVASVVRKEFDAQLQGNLKERNKFLAREASSLKYLAEEFSIPVILTNQITTHLSGALASQADLVSPADDLSLSEGTSGSSCVIAALGNTWSHSVNTRLILQYLDSERRQILIAKSPLAPFTSFVYTIKEEGLVLQAYGNS", 
#                 "MRGKTFRFEMQRDLVSFPLSPAVRVKLVSAGFQTAEELLEVKPSELSKEVGISKAEALETLQIIRRECLTNKPRYAGTSESHKKCTALELLEQEHTQGFIITFCSALDDILGGGVPLMKTTEICGAPGVGKTQLCMQLAVDVQIPECFGGVAGEAVFIDTEGSFMVDRVVDLATACIQHLQLIAEKHKGEEHRKALEDFTLDNILSHIYYFRCRDYTELLAQVYLLPDFLSEHSKVRLVIVDGIAFPFRHDLDDLSLRTRLLNGLAQQMISLANNHRLAVILTNQMTTKIDRNQALLVPALGESWGHAATIRLIFHWDRKQRLATLYKSPSQKECTVLFQIKPQGFRDTVVTSACSLQTEGSLSTRKRSRDPEEEL",
#                 "MGVLRVGLCPGLTEEMIQLLRSHRIKTVVDLVSADLEEVAQKCGLSYKALVALRRVLLAQFSAFPVNGADLYEELKTSTAILSTGIGSLDKLLDAGLYTGEVTEIVGGPGSGKTQVCLCMAANVAHGLQQNVLYVDSNGGLTASRLLQLLQAKTQDEEEQAEALRRIQVVHAFDIFQMLDVLQELRGTVAQQVTGSSGTVKVVVVDSVTAVVSPLLGGQQREGLALMMQLARELKTLARDLGMAVVVTNHITRDRDSGRLKPALGRSWSFVPSTRILLDTIEGAGASGGRRMACLAKSSRQPTGFQEMVDIGTWGTSEQSATLQGDQT",
#                 "GCSAFHRAESGTELLARLEGRSSLKEIEPNLFADEDSPVHGDILEFHGPEGTGKTEMLYHLTARCILPKSEGGLEVEVLFIDTDYHFDMLRLVTILEHRLSQSSEEIIKYCLGRFFLVYCSSSTHLLLTLYSLESMFCSHPSLCLLILDSLSAFYWIDRVNGGESVNLQESTLRKCSQCLEKLVNDYRLVLFATTQTIMQKASSSSEEPSHASRRLCDVDIDYRPYLCKAWQQLVKHRMFFSKQDDSQSSNQFSLVSRCLKSNSLKKHFFIIGESGVEFC"]
    
#     # 샘플 딕셔너리로 구성
#     samples = [
#         {"name": "T1109-dimer", "seqs": sequences1},
#         {"name": "H1166-trimer", "seqs": sequences2},
#         {"name": "H1185-tetramer", "seqs": sequences3},
#     ]

#     return render(request, 'nim/af2_multi/input.html', {
#         'databases': databases,
#         'samples_json': json.dumps(samples)
#     })

def nimAlphafoldMultiInputPage(request):
    databases = ["uniref90", "small_bfd", "mgnify"]

    return render(request, 'nim/af2_multi/input.html', {
        'databases': databases,
        'samples_json': json.dumps(multiSampleList())
    })

def nimAlphafoldMultiProcess(request):
    if request.method == "POST":    
        try:
            sample_dir = f"{outputSamplePath()}/nim/af2_multi"
            
            data = json.loads(request.body)
            sample_name = data.get('sampleName', '')    
            time.sleep(1)  # 5초 동안 멈춤
            results_dir = f"{sample_dir}/{sample_name}"
            return JsonResponse({"status": "done", "results_dir" : results_dir})  # ← 클라이언트에 저장 완료만 알림
        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)

def nimAlphafoldMultiListPage(request):
    try:
        results_dir = request.GET.get('results_dir')
        
        with open(f"{results_dir}/af2_multimer.json", 'r') as f:
            fold = json.load(f)
        if fold is None:
            return JsonResponse({"error": "세션에 저장된 결과가 없습니다."}, status=404)
        
        cards = []
        for i in range(len(fold)):
            title_pretty = f"Model {i + 1}"
            url = f"{reverse('nim:nimAlphafoldMultiOutputPage')}?results_dir={results_dir}&num={i}"
                
            cards.append({
                "title": title_pretty,
                "icon": "img/nvidia_alphafold2_multimer.jpg",
                "color": "#198754",
                "url": url,
                "delay": f"{0.1 + i * 0.1:.1f}s"
            })
        
        return render(request, 'list.html', {"cards": cards, "tool" : "Alphafold2 Multimer"})
    except Exception as e:
        print("에러 발생:", e)
        return JsonResponse({"error": str(e)}, status=500)  

def nimAlphafoldMultiOutputPage(request):
    results_dir = request.GET.get('results_dir')
    num = int(request.GET.get('num'))
    
    with open(f"{results_dir}/af2_multimer.json", 'r') as f:
        fold = json.load(f)
    if fold is None:
        return JsonResponse({"error": "세션에 저장된 결과가 없습니다."}, status=404)

    return render(request, 'nim/af2_multi/result.html', {"result_pdb" : fold[num], "num" : num + 1})







def nimRfdiffusionInputPage(request):
    # 샘플 딕셔너리로 구성 
    # samples = [
    #     {"pdbid": "1R42", "contigs": "A114-353/0 50-100", "hotspot_residues" : ["A119", "A123", "A233", "A234", "A235"]},
    #     {"pdbid": "5TPN", "contigs": "L1-25/0 70-100", "hotspot_residues" : ["L14", "L15", "L17", "L18"]},
    #     {"pdbid": "6VXX", "contigs": "A353-410/0 100-200", "hotspot_residues" : ["A360", "A361", "A362", "A366"]},
    # ]
    
    samples = [
        {"protein" : "TNIK", "pdbid": "2W3L", "contigs": "A9-32/0 49-165", "hotspot_residues" : ["A55", "A56", "A57", "A58", "A59"]},
        {"protein" : "BCL2", "pdbid": "2X7F", "contigs": "A15-60/0 100-313", "hotspot_residues" : ["A90", "A91", "A92", "A93"]},
        {"protein" : "HER2", "pdbid": "8VB5", "contigs": "A706-780/0 820-990", "hotspot_residues" : ["A800", "A801", "A802", "A803", "A804"]},
    ]

    return render(request, 'nim/rfd/input.html', {
        'samples_json': json.dumps(samples)
    })
    
def nimRfdiffusionProcess(request):
    if request.method == "POST":    
        try:
            sample_dir = f"{outputSamplePath()}/nim/rfd"
            
            data = json.loads(request.body)
            sample_name = data.get('sampleName', '')
            time.sleep(1)  # 5초 동안 멈춤     
            
            results_dir = f"{sample_dir}/{sample_name}"
            
            return JsonResponse({"status": "done", "results_dir" : results_dir})  # ← 클라이언트에 저장 완료만 알림
        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)
        
def nimRfdiffusionOutputPage(request):
    results_dir = request.GET.get('results_dir')
    
    # print("nimRfdiffusionOutputPage : ", sample_name)
    # rf디퓨션 : 디퓨션 스텝 50으로 돌린 결과를 불러옴
    sample_pdb_path = os.path.join(results_dir, "rfd.pdb")

    print(sample_pdb_path)
    
    pdb_content = ""
    if os.path.exists(sample_pdb_path):
        with open(sample_pdb_path, 'r') as f:
            pdb_content = f.read()

    print("pdb_content : ", pdb_content[:100])
    
    return render(request, 'nim/rfd/result.html', {
        "result_pdb": pdb_content
    })
    

def nimDiffdockInputPage(request):
    # 샘플 딕셔너리로 구성 
    samples = [
        {"protein" : "HIV-1 protease", "pdbid": "1HVR", "smiles": "CC(C)C[C@@H](NC(=O)[C@@H](Cc1ccccc1)NC(=O)[C@@H](Cc2ccccc2)NC(=O)C(C)(C)C)C(=O)N[C@@H](Cc3ccccc3)C(=O)N(C)C"},
        {"protein" : "Covid 19", "pdbid": "6LU7", "smiles": "CC(C)C[C@H](NC(=O)[C@H](CC1=CC=CC=C1)NC(=O)[C@H](CC2=CC=CC=C2)NC(=O)C(C)(C)C)C(=O)N[C@H](CC3=CC=CC=C3)C(=O)N(C)C"},
        # {"protein" : "EGFR kinase domain", "pdbid": "1M17", "smiles": "CN(C)CCOC1=NC2=C(C=C(C=C2)OC)N=C1C3=CC=CC=C3"},
        {"protein" : "Geldanamycin", "pdbid": "1YET", "smiles": "CC1CC(C2C(C1C(=O)NC3=CC=CC=C3)C(=O)C(=C(C2=O)OC)C)C(=O)C4=CC=CC=C4"},
        {"protein" : "Cyclooxygenase-2", "pdbid": "5IKR", "smiles": "CC1=C(C(=CC=C1)NC2=CC=CC=C2C(=O)O)C"},
    ]

    return render(request, 'nim/diffdock/input.html', {
        'samples_json': json.dumps(samples)
    })
    
def nimDiffdockProcess(request):
    if request.method == "POST":
        try:
            sample_dir = f"{outputPath()}/nim/diffdock/sample_result"
            data = json.loads(request.body)
            sample_name = data.get('sampleName', '')

            print(sample_name)
            
            results_dir = os.path.join(sample_dir, sample_name)
            time.sleep(1)  # 5초 동안 멈춤
            return JsonResponse({
                "status": "done",
                "results_dir": results_dir
            })

        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)

def nimDiffdockList(request):   
    try:
        results_dir = request.GET.get('results_dir')  # e.g., /path/to/sample_result/1HVR
        if not results_dir or not os.path.exists(results_dir):
            return JsonResponse({"error": "유효한 results_dir 경로가 아닙니다."}, status=400)

        pose_folders = sorted([
            f for f in os.listdir(results_dir)
            if os.path.isdir(os.path.join(results_dir, f)) and f.startswith("DiffdockPose")
        ])

        cards = []
        for i, folder_name in enumerate(pose_folders):
            title_pretty = f"Model {i + 1}"
            url = f"{reverse('nim:nimDiffdockOutputPage')}?results_dir={results_dir}/{folder_name}"

            cards.append({
                "title": title_pretty,
                "icon": "img/nvidia_diffdock.jpg",  # 추후 썸네일 이미지로 교체 가능
                "color": "#198754",
                "url": url,
                "delay": f"{0.1 + i * 0.1:.1f}s"
            })

        return render(request, 'list.html', {"cards": cards, "tool": "Diffdock"})
    
    except Exception as e:
        print("에러 발생:", e)
        return JsonResponse({"error": str(e)}, status=500)
    
def nimDiffdockOutputPage(request):
    results_dir = request.GET.get('results_dir')

    pdb_files = [f for f in os.listdir(results_dir) if f.endswith('.pdb')]
    
    pdb_path = os.path.join(results_dir, pdb_files[0])
    with open(pdb_path, 'r') as f:
        pdb_content = f.read()

    return render(request, 'nim/diffdock/result.html', {
        "result_pdb": pdb_content,  # JS에서 쓰게
        "results_dir": results_dir,
    })



def nimProteinmpnnInputPage(request):
    # 샘플 딕셔너리로 구성 
    # samples = [
    #     {"protein" : "HIV-1 protease", "pdbid": "1HVR"},
    #     {"protein" : "Covid 19", "pdbid": "6LU7"},
    #     {"protein" : "EGFR kinase domain", "pdbid": "1M17"},
    #     {"protein" : "Geldanamycin", "pdbid": "1YET"},
    #     {"protein" : "Cyclooxygenase-2", "pdbid": "5IKR"},
    # ]
    
    samples = [
        {"protein" : "TNIK", "pdbid": "2X7F"},
        {"protein" : "HER2", "pdbid": "8VB5"},
        {"protein" : "BCL2", "pdbid": "2W3L"},
    ]
    

    return render(request, 'nim/proteinmpnn/input.html', {
        'samples_json': json.dumps(samples)
    })
        
def nimProteinmpnnProcess(request):
    if request.method == "POST":    
        try:
            sample_dir = f"{outputSamplePath()}/nim/proteinmpnn_sample"
            
            data = json.loads(request.body)
            pdb_id = data.get('pdbId', '')
            temperature = data.get('temperature', '')
            
            print(pdb_id, temperature)
            results_dir = f"{sample_dir}/{pdb_id}/{pdb_id}_temp{''.join(temperature.split('.'))}"
            
            return JsonResponse({"status": "done", "results_dir" : results_dir})  # ← 클라이언트에 저장 완료만 알림
        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)

def nimProteinmpnnList(request):   
    try:
        results_dir = request.GET.get('results_dir')
        
        with open(f"{results_dir}/proteinmpnn_input.json", 'r') as f:
            fold = json.load(f)
        if fold is None:
            return JsonResponse({"error": "세션에 저장된 결과가 없습니다."}, status=404)
        
        cards = []
        for i in range(len(fold)):
            title_pretty = f"Model {i + 1}"
            url = f"{reverse('nim:nimProteinmpnnOutputPage')}?results_dir={results_dir}&num={i}"
                
            cards.append({
                "title": title_pretty,
                "icon": "img/nvidia_alphafold2_multimer.jpg",
                "color": "#198754",
                "url": url,
                "delay": f"{0.1 + i * 0.1:.1f}s"
            })

        return render(request, 'list.html', {"cards": cards, "tool": "ProteinMpnn"})
    except Exception as e:
        print("에러 발생:", e)
        return JsonResponse({"error": str(e)}, status=500)

def nimProteinmpnnOutputPage(request):
    results_dir = request.GET.get('results_dir')
    num = int(request.GET.get('num'))
    
    origin_pdb = f"{results_dir}/origin_protein.pdb"
    input_pdb_path = f"{results_dir}/proteinmpnn_input.json"
    pred_pdb_path = f"{results_dir}/proteinmpnn_temp.json"
    
    input_fa = f"{Path(results_dir).resolve().parent}/input.fa"
    temp_fa = f"{results_dir}/temp.fa"
    
    with open(origin_pdb, 'r') as f:
        origin_pdb_content = f.read()
    with open(input_pdb_path, 'r') as f:
        input_pdb_json = json.load(f)
    with open(pred_pdb_path, 'r') as f:
        temp_pdb_json = json.load(f)
    
    input_pdb_content = input_pdb_json[num][f'model{num+1}']
    print(input_pdb_content[:100])
    temp_pdb_content = temp_pdb_json[num][f'model{num+1}']
    print(temp_pdb_content[:100])
    
    with open(input_fa, 'r') as f:
        input_fa_content = f.read()
    with open(temp_fa, 'r') as f:
        temp_fa_content = f.read()
    
    # if openfold is None:
    #     return JsonResponse({"error": "세션에 저장된 결과가 없습니다."}, status=404)
    
    context = {
        "origin_pdb" : origin_pdb_content,
        "input_pdb" : input_pdb_content,
        "temp_pdb" : temp_pdb_content,
        "input_fa" : input_fa_content,
        "temp_fa" : temp_fa_content,
    }
    
    return render(request, "nim/proteinmpnn/result.html", context)