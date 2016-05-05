package interactome.analysis;


import interactome.Logger;
import interactome.data.BioDB;
import interactome.data.Refseq;
import interactome.input.Input;
import interactome.input.RefseqInput;

import org.rosuda.JRI.*;

public class BiasCorrector {
	static BiasCorrector _instance;
	
	Input input;
	DynamicParameters dynamicParameters;
	
	double alpha = 0;
	double beta = 0;
	
	public static BiasCorrector getInstance() {
		return _instance;
	}
	
	public static BiasCorrector createInstance(Input input, DynamicParameters dp) {
		_instance = new BiasCorrector(input, dp);
		return _instance;
	}

	private BiasCorrector(Input input, DynamicParameters dp) {
		this.input = input;
		this.dynamicParameters = dp;
	}
	
	public boolean correctBias() {
		BioDB biodb = BioDB.getInstance();
		
		Logger.logf("\ncorrecting bias.");
		
		// constants
		int max_frompolya = 3000;
		int iteration_count_max = 1000;
		double iteration_abort_ratio = 0.001;
		
		int n = this.dynamicParameters.parameterRefseqs.length;
		
		Refseq[] refseqs = new Refseq[n];
		RefseqInput[] inputs = new RefseqInput[n];
		double[] v = new double[n];
		double[] mappability_sums = new double[n];
		double[] weight = new double[n];
		double[] sums = new double[n];
		int[] mappable_position_counts = new int[n];
		int total_mappable_position_count = 0;
		
		// initialize
		for (int i=0; i<n; i++) {
			refseqs[i] = this.dynamicParameters.parameterRefseqs[i];
			inputs[i] = this.input.refseq_inputs.get(refseqs[i].refseq_id);
			
			for (int j=0; j<Math.min(max_frompolya, refseqs[i].length); j++) {
				mappability_sums[i] += refseqs[i].mappability[j];
				if (refseqs[i].mappability[j] > 0) {
					mappable_position_counts[i]++;
					total_mappable_position_count++;
				}
			}
			for (int j=0; j<Math.min(max_frompolya, refseqs[i].length); j++) {
				if (refseqs[i].mappability[j] > 0) {
					sums[i] += mappability_sums[i] / (mappable_position_counts[i] * refseqs[i].mappability[j]) * inputs[i].overlap_counts[j];
				} else if (inputs[i].overlap_counts[j] > 0) {
					Logger.errorf("unmappable position was overlapped by a read: %s %d", refseqs[i].refseq_id, j);
				}
			}
			
			v[i] = sums[i] / mappable_position_counts[i];
		}
		
		Rengine engine = new Rengine(new String[]{"--no-save"}, false, null);
		REXP result;
		
		// iteration
		for (int x=0; x<iteration_count_max; x++) {
			Logger.logf("iteration %d.", x+1);
			
			// solve alpha, beta by maximum-likelihood
			int[] X1 = new int[total_mappable_position_count];
			double[] X2 = new double[total_mappable_position_count];
			double[] X3 = new double[total_mappable_position_count];
			double[] offset = new double[total_mappable_position_count];
			
			// arrange data for R
			int count = 0;
			for (int i=0; i<n; i++) {
				for (int j=0; j<Math.min(max_frompolya, refseqs[i].length); j++) {
					if (refseqs[i].mappability[j] > 0) {
						X1[count] = (int)inputs[i].overlap_counts[j];
						X2[count] = ((double)refseqs[i].gc_percent[j]);
						X3[count] = ((double)refseqs[i].length - j);
						offset[count] = Math.log(v[i] * mappable_position_counts[i] * refseqs[i].mappability[j] / mappability_sums[i]);
						count++;
					}
				}
			}
			
			// execute R
			engine.assign("x1",X1);
			engine.assign("x2",X2);
			engine.assign("x3",X3);
			engine.assign("cs",offset);
			engine.eval("glmres<-glm(x1~x2+x3-1, family=poisson, offset=cs)");
			
			double old_alpha = alpha;
			double old_beta = beta;
			
			result = engine.eval("glmres$coefficients[1]");
			alpha = result.asDouble();
			result = engine.eval("glmres$coefficients[2]");
			beta = result.asDouble();
			Logger.logf("alpha: %.10f\nbeta: %.10f", alpha, beta);
			
			// change ratio
			boolean abort = false;
			if (Math.abs((alpha-old_alpha)/old_alpha) < iteration_abort_ratio
				&& Math.abs((beta - old_beta)/old_beta) < iteration_abort_ratio) {
				abort = true;
			}
			Logger.logf("change ratio -- alpha: %f, beta: %f, square: %f",
					Math.abs((alpha - old_alpha) / old_alpha),
					Math.abs((beta - old_beta) / old_beta),
					Math.sqrt(Math.pow(alpha-old_alpha, 2)+Math.pow(beta-old_beta, 2)) / Math.sqrt(alpha*alpha + beta*beta));
			
			// update v
			for (int i=0; i<n; i++) {
				weight[i] = 0;
				for (int j=0; j<Math.min(max_frompolya, refseqs[i].length); j++) {
					if (refseqs[i].mappability[j] > 0) {
						weight[i] += Math.exp(alpha * refseqs[i].gc_percent[j] + beta * (refseqs[i].length - j));
					}
				}
				v[i] = sums[i] / weight[i];
			}
			
			if (abort) break;
		}
		engine.end();
		
		// correction for all refseqs
		weight = new double[biodb.refseq_db.size()];
		mappable_position_counts = new int[biodb.refseq_db.size()];
		mappability_sums = new double[biodb.refseq_db.size()];
		sums = new double[biodb.refseq_db.size()];
		
		for (int i=0; i<biodb.all_refseq_ids.length; i++) {
			String refseq_id = biodb.all_refseq_ids[i];
			Refseq refseq = biodb.refseq_db.get(refseq_id);
			RefseqInput refinput = input.refseq_inputs.get(refseq_id);
			
			if (refseq != null && !refseq.is_invalid) {
				double w = 0;
				for (int j=0; j<refseq.length; j++) {
					if (refseq.mappability[j] > 0) {
						w += Math.exp(alpha * refseq.gc_percent[j] + beta * (refseq.length - j));
					}
				}
				
				for (int j=0; j<refseq.length; j++) {
					mappability_sums[i] += refseq.mappability[j];
					if (refseq.mappability[j] > 0) {
						mappable_position_counts[i]++;
					}
				}
				
				for (int j=0; j<refseq.length; j++) {
					if (refseq.mappability[j] > 0) {
						sums[i] += mappability_sums[i] / (mappable_position_counts[i] * refseq.mappability[j]) * refinput.starting_counts[j];
					}
				}
				
				refinput.true_expression = sums[i] / w;
			}
		}
		Logger.logf("bias correction done.");
		
		return true;
	}
}
