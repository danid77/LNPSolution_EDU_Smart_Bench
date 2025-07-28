import os, sys
from datetime import datetime
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

from property_cal_founction import *

if __name__ == "__main__":    
    # current_time = datetime.now().strftime("%Y%m%d_%Hh%Mm%Ss")
    
    if len(sys.argv) < 2:
        print("Usage: <csv files 이름> <출력할 umap 갯수>")
        sys.exit(1)
        
    in_csv_file = sys.argv[1].strip('"').strip()
    out_csv_file_path = sys.argv[2].strip('"').strip()
    cnt = 2
    
    current_dir = os.getcwd()
    csv_file = os.path.join(current_dir, in_csv_file)
    csv_file_name = in_csv_file.split(".")[0]
    
    input_df = pd.read_csv(csv_file)
    input_df.loc[: , "Status"] = "Generate"

    # 선택할 컬럼 리스트 (우선순위: SMILES > Smiles > smiles)
    column_options = ["SMILES", "Smiles", "smiles"]

    # 존재하는 컬럼 중 하나 선택
    selected_column = next((col for col in column_options if col in input_df.columns), None)
    # print(input_df)
    # print(input_df.columns)
    
    
    if selected_column:
        _input_smiles = input_df.at[0, 'Input_SMILES']
        input_row = pd.DataFrame([{selected_column : _input_smiles, 'Tanimoto' : 1, 'NLL' : 0, 'Status' : "Query"}])
        input_df = input_df.drop(columns=['Input_SMILES'])
        # concat(axis=0): 행 기준으로 합치기
        df = pd.concat([input_row, input_df], ignore_index=True)
        # print(df)
        smiles = df[[selected_column]].dropna()
    else:
        raise ValueError("SMILES, Smiles, smiles 컬럼이 존재하지 않습니다.")
    
    # 🔹 0~2047 컬럼 중 하나라도 NaN이면 삭제
    fingerprint_columns = [f"fp_{i+1}" for i in range(1, 2048)]  # 0~2047까지 컬럼 리스트  
    _cal_result = process_smiles_file_ecfp6(smiles, selected_column)
    cal_result = _cal_result.dropna(subset=fingerprint_columns)
    before_drop, after_drop = len(_cal_result), len(cal_result)
    print(cal_result.head())
    print(f"\n🔹 총 {before_drop}개 중 {after_drop}개가 남음. ({before_drop - after_drop} 삭제 됨)")
    
    # 좌표 만들기
    # umap_df = umap_coordinate(cal_result, selected_column, cnt)
    umap_df = umap_coordinate(_cal_result, selected_column, cnt)
    smiles_umap_df = pd.concat([df, umap_df], axis=1)

    # tsne_df = tsne_coordinate(cal_result, selected_column, cnt)
    # smiles_tsne_df = pd.concat([df, tsne_df], axis=1)

    cal_df = pd.concat([df.drop(columns=selected_column), cal_result], axis=1)
    print(cal_df.head())
    cal_df.to_csv(f"{out_csv_file_path}/Rdkit_calculate_result.csv", index=False)
    
    smiles_umap_df.to_csv(f"{out_csv_file_path}/umap_result.csv", index=False)
    # smiles_tsne_df.to_csv(f"{out_csv_file_path}/tSNE_result.csv", index=False)
    
    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        x=smiles_umap_df['UMAP1'], 
        y=smiles_umap_df['UMAP2'], 
        hue=smiles_umap_df['Status'], 
        palette='tab10', 
        alpha=0.7, 
        s=50
    )
    plt.xlabel("UMAP 1")
    plt.ylabel("UMAP 2")
    plt.title("UMAP Scatter Plot (Colored by data_type)")
    plt.legend(title="Status", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True)

    # ✅ 이미지 파일로 저장
    plt.savefig(f"{out_csv_file_path}/umap.png", dpi=300, bbox_inches='tight')
    
    # plt.figure(figsize=(8, 6))
    # sns.scatterplot(x=smiles_tsne_df['tSNE1'], y=smiles_tsne_df['tSNE2'], hue=smiles_tsne_df['tool'], palette='tab10', alpha=0.7, s=50)
    # plt.xlabel("tSNE 1")
    # plt.ylabel("tSNE 2")
    # plt.title("tSNE Scatter Plot (Colored by tool)")
    # plt.legend(title="tool", bbox_to_anchor=(1.05, 1), loc='upper left')
    # plt.grid(True)
    
    # # ✅ 이미지 파일로 저장
    # plt.savefig(f"{out_csv_file_path}/tsne.png", dpi=300, bbox_inches='tight')