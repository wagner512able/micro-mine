# benchmarks/scaling.py
import time
import tracemalloc
import random
import re
import os

def generate_synthetic_genome(size_bp=155055, substitution_rate=0.005):
    """Gera um genoma sintético de ~155kb com substituições aleatórias."""
    bases = ['A', 'T', 'C', 'G']
    # Sequência base simples (pode ser aleatória)
    seq = [random.choice(bases) for _ in range(size_bp)]
    
    # Aplica substituições para simular divergência entre genomas
    num_subs = int(size_bp * substitution_rate)
    for _ in range(num_subs):
        idx = random.randint(0, size_bp - 1)
        seq[idx] = random.choice(bases)
        
    return "".join(seq)

def mine_ssr_simple(sequence):
    """Simula a mineração de SSRs (mono a hexa) usando regex."""
    # Padrões simplificados para benchmark (MISA defaults)
    patterns = [
        r'(A{10,}|T{10,}|C{10,}|G{10,})', # Mono
        r'([ATCG]{2}){6,}', # Di
        r'([ATCG]{3}){5,}', # Tri
        r'([ATCG]{4}){5,}', # Tetra
        r'([ATCG]{5}){5,}', # Penta
        r'([ATCG]{6}){5,}'  # Hexa
    ]
    ssrs_found = 0
    for p in patterns:
        ssrs_found += len(re.findall(p, sequence))
    return ssrs_found

def run_scaling_benchmark():
    print("Iniciando benchmark de escalabilidade (benchmarks/scaling.py)...")
    print(f"{'Genomas':<10} | {'SSRs':<10} | {'Tempo (s)':<10} | {'Memória Pico (MB)':<15}")
    print("-" * 55)
    
    genome_sizes = [2, 5, 10, 25, 50]
    
    for n_genomes in genome_sizes:
        # Inicia medição de memória
        tracemalloc.start()
        start_time = time.time()
        
        total_ssrs = 0
        genomes = []
        
        # Gera N genomas sintéticos
        for i in range(n_genomes):
            # 0.5% de substituição para cada cópia, como no manuscrito
            seq = generate_synthetic_genome(size_bp=155055, substitution_rate=0.005)
            genomes.append(seq)
            
        # Roda a mineração em todos os genomas
        for seq in genomes:
            total_ssrs += mine_ssr_simple(seq)
            
        end_time = time.time()
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        
        peak_mb = peak / (1024 * 1024)
        elapsed = end_time - start_time
        
        print(f"{n_genomes:<10} | {total_ssrs:<10} | {elapsed:<10.2f} | {peak_mb:<15.2f}")

if __name__ == "__main__":
    run_scaling_benchmark()