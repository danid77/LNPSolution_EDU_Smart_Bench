from django.shortcuts import render
from django.http import JsonResponse
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

from apps import getKey, outputPath, outputSamplePath, sampleList
from apps.umap.services import molMin, genMol, molsparkProcess
# Create your views here.
        
def nimOpenfoldInputPageSample(request):
    databases = ['Uniref30_2302', 'PDB70_220313', 'colabfold_envdb_202108']
    
    sequences1 = "MATKAVCVLKGDGPVQGIINFEQKESNGPVKVWGSIKGLTEGLHGFHVHEFGDNTAGCTSAGPHFNPLSRKHGGPKDEERHVGDLGNVTADKDGVADVSIEDSVISLSGDHCIIGRTLVVHEKADDLGKGGNEESTKTGNAGSRLACGVIGIAQMATKAVCVLKGDGPVQGIINFEQKESNGPVKVWGSIKGLTEGLHGFHVHEFGDNTAGCTSAGPHFNPLSRKHGGPKDEERHVGDLGNVTADKDGVADVSIEDSVISLSGDHCIIGRTLVVHEKADDLGKGGNEESTKTGNAGSRLACGVIGIAQ"
    sequences2 = "MAHAAQVGLQDATSPIMEELITFHDHALMIIFLICFLVLYALFLTLTTKLTNTNISDAQEMETVWTILPAIILVLIALPSLRILYMTDEVNDPSLTIKSIGHQWYWTYEYTDYGGLIFNSYMLPPLFLEPGDLRLLDVDNRVVLPIEAPIRMMITSQDVLHSWAVPTLGLKTDAIPGRLNQTTFTATRPGVYYGQCSEICGANHSFMPIVLELIPLKIFEMGPVFTL"
    sequences3 = "MSTESMIRDVELAEEALPKKTGGPQGSRRCLFLSLFSFLIVAGATTLFCLLHFGVIGPQREEFPRDLSLISPLAQAVRSSSRTPSDKPVAHVVANPQAEGQLQWLNRRANALLANGVELRDNQLVVPSEGLYLIYSQVLFKGQGCPSTHVLLTHTISRIAVSYQTKVNLLSAIKSPCQRETPEGAEAKPWYEPIYLGGVFQLEKGDRLSAEINRPDYLDFAESGQVYFGIIALMSTESMIRDVELAEEALPKKTGGPQGSRRCLFLSLFSFLIVAGATTLFCLLHFGVIGPQREEFPRDLSLISPLAQAVRSSSRTPSDKPVAHVVANPQAEGQLQWLNRRANALLANGVELRDNQLVVPSEGLYLIYSQVLFKGQGCPSTHVLLTHTISRIAVSYQTKVNLLSAIKSPCQRETPEGAEAKPWYEPIYLGGVFQLEKGDRLSAEINRPDYLDFAESGQVYFGIIALMSTESMIRDVELAEEALPKKTGGPQGSRRCLFLSLFSFLIVAGATTLFCLLHFGVIGPQREEFPRDLSLISPLAQAVRSSSRTPSDKPVAHVVANPQAEGQLQWLNRRANALLANGVELRDNQLVVPSEGLYLIYSQVLFKGQGCPSTHVLLTHTISRIAVSYQTKVNLLSAIKSPCQRETPEGAEAKPWYEPIYLGGVFQLEKGDRLSAEINRPDYLDFAESGQVYFGIIAL"
    sequences4 = "MADSRDPASDQMQHWKEQRAAQKADVLTTGAGNPVGDKLNVITVGPRGPLLVQDVVFTDEMAHFDRERIPERVVHAKGAGAFGYFEVTHDITKYSKAKVFEHIGKKTPIAVRFSTVAGESGSADTVRDPRGFAVKFYTEDGNWDLVGNNTPIFFIRDPILFPSFIHSQKRNPQTHLKDPDMVWDFWSLRPESLHQVSFLFSDRGIPDGHRHMNGYGSHTFKLVNANGEAVYCKFHYKTDQGIKNLSVEDAARLSQEDPDYGIRDLFNAIATGKYPSWTFYIQVMTFNQAETFPFNPFDLTKVWPHKDYPLIPVGKLVLNRNPVNYFAEVEQIAFDPSNMPPGIEASPDKMLQGRLFAYPDTHRHRLGPNYLHIPVNCPYRARVANYQRDGPMCMQDNQGGAPNYYPNSFGAPEQQPSALEHSIQYSGEVRRFNTANDDNVTQVRAFYVNVLNEEQRKRLCENIAGHLKDAQIFIQKKAVKNFTEVHPDYGSHIQALLDKYNAEKPKNAIHTFVQSGSHLAAREKANL"
    sequences5 = "MDGETAEEQGGPVPPPVAPGGPGLGGAPGGRREPKKYAVTDDYQLSKQVLGLGVNGKVLECFHRRTGQKCALKLLYDSPKARQEVDHHWQASGGPHIVCILDVYENMHHGKRCLLIIMECMEGGELFSRIQERGDQAFTEREAAEIMRDIGTAIQFLHSHNIAHRDVKPENLLYTSKEKDAVLKLTDFGFAKETTQNALQTPCYTPYYVAPEVLGPEKYDKSCDMWSLGVIMYILLCGFPPFYSNTGQAISPGMKRRIRLGQYGFPNPEWSEVSEDAKQLIRLLLKTDPTERLTITQFMNHPWINQSMVVPQTPLHTARVLQEDKDHWDEVKEEMTSALATMRVDYDQVKIKDLKTSNNRLLNKRRKKQAGSSSASQGCNNQ"
    
    #     # 샘플 딕셔너리로 구성 
    # samples = [
    #     {"protein" : "HIV-1 protease", "pdbid": "1HVR"},
    #     {"protein" : "Covid 19", "pdbid": "6LU7"},
    #     {"protein" : "EGFR kinase domain", "pdbid": "1M17"},
    #     {"protein" : "Geldanamycin", "pdbid": "1YET"},
    #     {"protein" : "Cyclooxygenase-2", "pdbid": "5IKR"},
    # ]
    
    # 샘플 딕셔너리로 구성
    samples = [
        {"name": "SODC", "seq": sequences1},
        {"name": "COX2", "seq": sequences2},
        {"name": "TNFA", "seq": sequences3},
        {"name": "CATA", "seq": sequences4},
        {"name": "MAPK3", "seq": sequences5},
    ]
    
    return render(request, 'nim/openfold/input_sample.html', {'databases': databases, 'samples_json': json.dumps(samples)})

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

def nimOpenfoldOutputPage(request):
    results_dir = request.GET.get('results_dir')
    with open(f"{results_dir}/openfold.json", 'r') as f:
        openfold = json.load(f)
    if openfold is None:
        return JsonResponse({"error": "세션에 저장된 결과가 없습니다."}, status=404)
    
    result_pdb = openfold['structures_in_ranked_order'][0]['structure']
    confidence = openfold['structures_in_ranked_order'][0]['confidence']
    return render(request, 'nim/openfold/result.html', {'result_pdb' : result_pdb, "confidence" : confidence})


def nimMsaInputPage(request):
    sections = ["OpenFold2", "MMseqs2 MSA", "MolMin", "DiffDock Inputs"]
    databases = ["Uniref30_2303", "PDB70_230113", "colabfold_envdb_202108"]
    return render(request, 'nim/msa/input.html', {
        'sections': sections,
        'databases': databases
    })



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




def nimAlphafoldInputPage(request):
    databases = ["uniref90", "small_bfd", "mgnify"]
    df = pd.read_csv(f"{outputSamplePath()}/nim/af2/af2_sample_up.csv")
    sample_seq = df['seqeunce'].tolist()
    sample_name = df['name'].tolist()
        # zip으로 묶어서 전달
    samples = list(zip(sample_seq, sample_name))

    return render(request, 'nim/af2/input.html', {
        'databases': databases,
        'samples': samples  # <- 이거만 넘기면 됨
    })
    
def nimAlphafoldProcess(request):
    if request.method == "POST":    
        try:
            results_dir = f"{outputSamplePath()}/nim/af2"
        
            data = json.loads(request.body)
            sample_name = data.get('sampleName', '')   
            time.sleep(1)  # 5초 동안 멈춤
            return JsonResponse({"status": "done", "results_dir" : results_dir, "sample_name" : sample_name})  # ← 클라이언트에 저장 완료만 알림
        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)
        
def nimAlphafoldOutputPage(request):
    sample_name = request.GET.get('sample_name')
    results_dir = request.GET.get('results_dir')
    
    af2_result = ""
    if sample_name == "Sample1 PTEN":
        with open(f'{results_dir}/sample1_pten.json', 'r') as f1:
            af2_result = json.load(f1)
    elif sample_name == "Sample2 KRAS":
        with open(f'{results_dir}/sample2_kras.json', 'r') as f2:
            af2_result = json.load(f2)
    elif sample_name == "Sample3 CDK4":
        with open(f'{results_dir}/sample3_cdk.json', 'r') as f3:
            af2_result = json.load(f3) 
    return render(request, 'nim/af2/result.html', {"result_pdbs" : af2_result})

def nimAlphafoldMultiInputPage(request):
    databases = ["uniref90", "small_bfd", "mgnify"]
    
    sequences1 = ["MGSKKLKRVGLSQELCDRLSRHQILTCQDFLCLSPLELMKVTGLSYRGVHELLCMVSRACAPKMQTAYGIKAQRSADFSPAFLSTTLSALDEALHGGVACGSLTEITGPPGCGKTQFCIMMSILATLPTNMGGLEGAVVYIDTESAFSAERLVEIAESRFPRYFNTEEKLLLTSSKVHLYRELTCDEVLQRIESLEEEIISKGIKLVILDSVASVVRKEFDAQLQGNLKERNKFLAREASSLKYLAEEFSIPVILTNQITTHLSGALASQADLVSPADDLSLSEGTSGSSCVIAALGNTWSHSVNTRLILQYLDSERRQILIAKSPLAPFTSFVYTIKEEGLVLQAYGNS", 
              "GSHMAQPRPPFHITIPIYPGVDLLDVAAPVELFSWMADAWKARATTITLAAEHLTPLKTRDGLTLTPQRQFADYADAAAPQPQTHLLWVPGGAPDVLRKLMRGGPYLDFLKAQSAGADHVSSVCEGALLLAAAGLLDGYRATTHWAFIPCLQQFPAIKVAEGFPRYVIDGNRITGGGISSGLAEALAIVARVAGQDIAKHVQMITQYFPDPPFEQTIVPATHCPLQA"]  
    sequences2 = ["QSALTQPASVSGSPGQSITISCTGTSSDVGGYNYVSWYQQHPGKAPKLMIYEVTNRPSGVSNRFSGSRSGNTASLTISGLQAEDEADYYCSSYTSSSLYVFGTGTKVAVLGQPKANPTVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADSSPVKAGVETTKPSKQSNNKYAASSYLSLTPEQWKSHRSYSCQVTHEGSTVEKTVAPTECS", 
                "QVHLVQSGAEVKKPGSSVKVSCKASGGTFSSCAISWVRQAPGQGLEWMGRIIPILGIANYAQKFQGRVTITADKSTSTAYMELSSLRSEDTAVYYCARGWEFGSGSYYRTDYYYYAMDVWGQGTTVTVSSASTKGPSVFPLAPCSRSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPK",
                "SNMNNTASWFTALTQHGKEDLKFPRGQGVPINTNSSPDDQIGYYRRATRRIRGGDGKMKDLSPRWYFYYLGTGPEAGLPYGANKDGIIWVATEGALNTPKDHIGTRNPANNAAIVLQLPQGTTLPKGFYA"]  
    sequences3 = ["MGSKKLKRVGLSQELCDRLSRHQILTCQDFLCLSPLELMKVTGLSYRGVHELLCMVSRACAPKMQTAYGIKAQRSADFSPAFLSTTLSALDEALHGGVACGSLTEITGPPGCGKTQFCIMMSILATLPTNMGGLEGAVVYIDTESAFSAERLVEIAESRFPRYFNTEEKLLLTSSKVHLYRELTCDEVLQRIESLEEEIISKGIKLVILDSVASVVRKEFDAQLQGNLKERNKFLAREASSLKYLAEEFSIPVILTNQITTHLSGALASQADLVSPADDLSLSEGTSGSSCVIAALGNTWSHSVNTRLILQYLDSERRQILIAKSPLAPFTSFVYTIKEEGLVLQAYGNS", 
                "MRGKTFRFEMQRDLVSFPLSPAVRVKLVSAGFQTAEELLEVKPSELSKEVGISKAEALETLQIIRRECLTNKPRYAGTSESHKKCTALELLEQEHTQGFIITFCSALDDILGGGVPLMKTTEICGAPGVGKTQLCMQLAVDVQIPECFGGVAGEAVFIDTEGSFMVDRVVDLATACIQHLQLIAEKHKGEEHRKALEDFTLDNILSHIYYFRCRDYTELLAQVYLLPDFLSEHSKVRLVIVDGIAFPFRHDLDDLSLRTRLLNGLAQQMISLANNHRLAVILTNQMTTKIDRNQALLVPALGESWGHAATIRLIFHWDRKQRLATLYKSPSQKECTVLFQIKPQGFRDTVVTSACSLQTEGSLSTRKRSRDPEEEL",
                "MGVLRVGLCPGLTEEMIQLLRSHRIKTVVDLVSADLEEVAQKCGLSYKALVALRRVLLAQFSAFPVNGADLYEELKTSTAILSTGIGSLDKLLDAGLYTGEVTEIVGGPGSGKTQVCLCMAANVAHGLQQNVLYVDSNGGLTASRLLQLLQAKTQDEEEQAEALRRIQVVHAFDIFQMLDVLQELRGTVAQQVTGSSGTVKVVVVDSVTAVVSPLLGGQQREGLALMMQLARELKTLARDLGMAVVVTNHITRDRDSGRLKPALGRSWSFVPSTRILLDTIEGAGASGGRRMACLAKSSRQPTGFQEMVDIGTWGTSEQSATLQGDQT",
                "GCSAFHRAESGTELLARLEGRSSLKEIEPNLFADEDSPVHGDILEFHGPEGTGKTEMLYHLTARCILPKSEGGLEVEVLFIDTDYHFDMLRLVTILEHRLSQSSEEIIKYCLGRFFLVYCSSSTHLLLTLYSLESMFCSHPSLCLLILDSLSAFYWIDRVNGGESVNLQESTLRKCSQCLEKLVNDYRLVLFATTQTIMQKASSSSEEPSHASRRLCDVDIDYRPYLCKAWQQLVKHRMFFSKQDDSQSSNQFSLVSRCLKSNSLKKHFFIIGESGVEFC"]
    
    # 샘플 딕셔너리로 구성
    samples = [
        {"name": "T1109-dimer", "seqs": sequences1},
        {"name": "H1166-trimer", "seqs": sequences2},
        {"name": "H1185-tetramer", "seqs": sequences3},
    ]

    return render(request, 'nim/af2_multi/input.html', {
        'databases': databases,
        'samples_json': json.dumps(samples)
    })
    
def nimAlphafoldMultiProcess(request):
    if request.method == "POST":    
        try:
            results_dir = f"{outputSamplePath()}/nim/af2_multi"
            
            data = json.loads(request.body)
            sample_name = data.get('sampleName', '')    
            time.sleep(1)  # 5초 동안 멈춤
            return JsonResponse({"status": "done", "results_dir" : results_dir, "sample_name" : sample_name})  # ← 클라이언트에 저장 완료만 알림
        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)
        
def nimAlphafoldMultiOutputPage(request):
    sample_name = request.GET.get('sample_name')
    results_dir = request.GET.get('results_dir')
    
    print(sample_name)
    sample_pdb_path = os.path.join(results_dir, f"{sample_name}.pdb")

    pdb_content = ""
    if os.path.exists(sample_pdb_path):
        with open(sample_pdb_path, 'r') as f:
            pdb_content = f.read()

    return render(request, 'nim/af2_multi/result.html', {
        "result_pdbs": json.dumps([pdb_content])
    })

  
def nimEsmfoldInputPage(request):
    df = pd.read_csv(f"{outputPath()}/nim/af2/af2_sample_up.csv")
    sample_seq = df['seqeunce'].tolist()
    sample_name = df['name'].tolist()
        # zip으로 묶어서 전달
    samples = list(zip(sample_seq, sample_name))

    return render(request, 'nim/esmfold/input.html', {
        'samples': samples  # <- 이거만 넘기면 됨
    })

def nimEsmfoldProcess(request):
    if request.method == "POST":    
        try:
            results_dir = f"{outputSamplePath()}/nim/esmfold"
        
            data = json.loads(request.body)
            sample_name = data.get('sampleName', '')   
            time.sleep(1)  # 5초 동안 멈춤
            return JsonResponse({"status": "done", "results_dir" : results_dir, "sample_name" : sample_name})  # ← 클라이언트에 저장 완료만 알림
        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)

def nimEsmOutputPage(request):
    sample_name = request.GET.get('sample_name')
    results_dir = request.GET.get('results_dir')
    
    esm_result = ""
    if sample_name == "Sample1 PTEN":
        with open(f'{results_dir}/sample1_pten.json', 'r') as f1:
            esm_result = json.load(f1)
    elif sample_name == "Sample2 KRAS":
        with open(f'{results_dir}/sample2_kras.json', 'r') as f2:
            esm_result = json.load(f2)
    elif sample_name == "Sample3 CDK4":
        with open(f'{results_dir}/sample3_cdk.json', 'r') as f3:
            esm_result = json.load(f3) 
            
    return render(request, 'nim/esmfold/result.html', {'result_pdb' : esm_result["pdbs"]})


def nimRfdiffusionInputPage(request):
    # 샘플 딕셔너리로 구성 
    samples = [
        {"pdbid": "1R42", "contigs": "A114-353/0 50-100", "hotspot_residues" : ["A119", "A123", "A233", "A234", "A235"]},
        {"pdbid": "5TPN", "contigs": "L1-25/0 70-100", "hotspot_residues" : ["L14", "L15", "L17", "L18"]},
        {"pdbid": "6VXX", "contigs": "A353-410/0 100-200", "hotspot_residues" : ["A360", "A361", "A362", "A366"]},
    ]

    return render(request, 'nim/rfd/input.html', {
        'samples_json': json.dumps(samples)
    })
    
def nimRfdiffusionProcess(request):
    if request.method == "POST":    
        try:
            results_dir = f"{outputPath()}/nim/rfd"
            
            data = json.loads(request.body)
            sample_name = data.get('sampleName', '')
            time.sleep(1)  # 5초 동안 멈춤     
            
            return JsonResponse({"status": "done", "results_dir" : results_dir, "sample_name" : sample_name})  # ← 클라이언트에 저장 완료만 알림
        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)
        
def nimRfdiffusionOutputPage(request):
    sample_name = request.GET.get('sample_name')
    results_dir = request.GET.get('results_dir')
    
    print("nimRfdiffusionOutputPage : ", sample_name)
    # rf디퓨션 : 디퓨션 스텝 50으로 돌린 결과를 불러옴
    sample_pdb_path = os.path.join(results_dir, f"rf_{sample_name}.pdb")

    pdb_content = ""
    if os.path.exists(sample_pdb_path):
        with open(sample_pdb_path, 'r') as f:
            pdb_content = f.read()

    return render(request, 'nim/rfd/result.html', {
        "result_pdbs": json.dumps([pdb_content])
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
            results_dir = f"{outputPath()}/nim/diffdock/sample_result"
            data = json.loads(request.body)
            sample_name = data.get('sampleName', '')

            sample_pdb_path = os.path.join(results_dir, sample_name)
            sample_list = sorted([
                name for name in os.listdir(sample_pdb_path)
                if os.path.isdir(os.path.join(sample_pdb_path, name))
            ])
            time.sleep(1)  # 5초 동안 멈춤
            return JsonResponse({
                "status": "done",
                "results_dir": sample_pdb_path,
                "sample_list": sample_list
            })

        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)
        
def nimDiffdockOutputPage(request):
    results_dir = request.GET.get('results_dir')
    sample_list_raw = request.GET.get('sample_list')  # JSON string 형태

    # URL 디코딩 먼저
    results_dir = unquote(results_dir)

    try:
        sample_list = json.loads(sample_list_raw)
    except Exception as e:
        sample_list = []

    return render(request, 'nim/diffdock/result.html', {
        "sample_list": json.dumps(sample_list),  # JS에서 쓰게
        "results_dir": results_dir,
    })



def nimProteinmpnnInputPage(request):
    # 샘플 딕셔너리로 구성 
    samples = [
        {"protein" : "HIV-1 protease", "pdbid": "1HVR"},
        {"protein" : "Covid 19", "pdbid": "6LU7"},
        {"protein" : "EGFR kinase domain", "pdbid": "1M17"},
        {"protein" : "Geldanamycin", "pdbid": "1YET"},
        {"protein" : "Cyclooxygenase-2", "pdbid": "5IKR"},
    ]

    return render(request, 'nim/proteinmpnn/input.html', {
        'samples_json': json.dumps(samples)
    })
        
def nimProteinmpnnProcess(request):
    if request.method == "POST":    
        try:
            sample_dir = f"{outputPath()}/nim/proteinmpnn_sample"
            
            data = json.loads(request.body)
            pdb_id = data.get('pdbId', '')
            print(pdb_id)
            results_dir = f"{sample_dir}/{pdb_id}"
            
            return JsonResponse({"status": "done", "results_dir" : results_dir})  # ← 클라이언트에 저장 완료만 알림
        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)
        
def nimProteinmpnnOutputPage(request):
    results_dir = request.GET.get('results_dir')
    pdb_id = os.path.basename(results_dir).lower()
    input_pdb_path = f"{results_dir}/{pdb_id}_input_model.pdb"
    pred_pdb_path = f"{results_dir}/{pdb_id}_t_model.pdb"
    mfasta_path = f"{results_dir}/{pdb_id}.fasta"
    
    with open(input_pdb_path, "r") as f1:
        input_pdb = f1.read()
    with open(pred_pdb_path, "r") as f2:
        pred_pdb = f2.read()
    with open(mfasta_path, "r") as f3:
        mfasta = f3.read()
        
    if input_pdb_path is None:
        return JsonResponse({"error": "세션에 저장된 결과가 없습니다."}, status=404)
    
    return render(request, "nim/proteinmpnn/result.html", {
        "input_pdb" : input_pdb, "pred_pdb" : pred_pdb, "result_mfasta" : mfasta
    })