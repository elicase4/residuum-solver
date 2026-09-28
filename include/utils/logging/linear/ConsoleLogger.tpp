namespace residuum::utils::logging::linear {

	template<typename DataT>
	void ConsoleLogger::log(Index iter, const std::vector<DataT>& perDOFAbs, DataT flopsThisIter) const {

		// capture the relative-residual baseline and start the clock on the first call
		if (!baselineCaptured_) {
			initialPerDOFNorms_.assign(perDOFAbs.begin(), perDOFAbs.end());
			baselineCaptured_ = true;
			startTime_ = std::chrono::steady_clock::now();
		}

		// lifetime accounting, used only by summary()
		totalFlops_ += static_cast<double>(flopsThisIter);
		lastIter_ = iter;

		// convert to relative, per-DOF, against the iteration-0 baseline
		lastRelPerDOF_.assign(perDOFAbs.size(), Real(0));
		for (Index i = 0; i < perDOFAbs.size(); ++i) {
			Real base = initialPerDOFNorms_[i];
			lastRelPerDOF_[i] = (base > Real(0)) ? static_cast<Real>(perDOFAbs[i]) / base : Real(0);
		}

		const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime_).count();

		// CSV row: every iteration, independent of the console interval throttle below
		if (csv_->enabled()) {
			std::vector<Real> row = {static_cast<Real>(outerTick_), static_cast<Real>(iter)};
			for (const auto& v : perDOFAbs) row.push_back(static_cast<Real>(v));
			for (const auto& v : lastRelPerDOF_) row.push_back(v);
			row.push_back(static_cast<Real>(flopsThisIter));
			row.push_back(static_cast<Real>(elapsed));
			csv_->writeRow(row);
		}

		// check to print based on config
		if (interval == 0) return;
		if ((iter > 0) && ((iter % interval) != 0)) return;

		// print header
		if (printHeader_ && iter == 0) {
			printBanner();
			printColumnHeader();
			printHeader_ = false;
		}

		// set print format
		*out_ << tag() << std::left << std::setw(6) << iter << "  ";
		*out_ << std::scientific << std::setprecision(4);

		for (const auto& v : lastRelPerDOF_) {
			*out_ << std::setw(14) << v << "  ";
		}

		*out_ << std::fixed << std::setprecision(3) << elapsed << "s";
		*out_ << "\n";

	}

	inline void ConsoleLogger::summary(bool converged) const {

		const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime_).count();
		const double avgFlopsPerIter = (lastIter_ > 0) ? (totalFlops_ / static_cast<double>(lastIter_)) : 0.0;

		*out_ << tag() << (converged ? "converged" : "did not converge") << ": " << lastIter_ << " iters, ";

		for (Index i = 0; i < lastRelPerDOF_.size(); ++i) {
			*out_ << "res[" << dofNames[i] << "] " << std::scientific << std::setprecision(4) << lastRelPerDOF_[i];
			if (i + 1 < lastRelPerDOF_.size()) *out_ << ", ";
		}

		*out_ << ", avg " << std::scientific << std::setprecision(2) << avgFlopsPerIter << " flops/iter, " << std::fixed << std::setprecision(3) << elapsed << "s\n";

	}

	inline void ConsoleLogger::printBanner() const {

		const int width = 60;

		*out_ << "\n";
		*out_ << tag() << std::string(width, '=') << "\n";
		*out_ << tag() << "(" << preconditionerName << ")\n";

		for (const auto& param : extraParams) {
			*out_ << tag() << param.first << ": " << param.second << "\n";
		}

		if (dofNames.size() > 1) {
			*out_ << tag() << "Components:";
			for (const auto& n : dofNames) *out_ << "  " << n;
			*out_ << "\n";
			*out_ << tag() << "DOF ordering: " << (dofOrdering == fem::dof::DOFOrdering::Interleaved ? "Interleaved (node-major)" : "Block (field-major)") << "\n";
		}

		*out_ << tag() << std::string(width, '=') << "\n";

	}

	inline void ConsoleLogger::printColumnHeader() const {

		*out_ << tag() << std::left << std::setw(6) << "Iter" << "  ";

		for (const auto& name : dofNames) {
			*out_ << std::setw(14) << ("Res[" + name + "]") << "  ";
		}

		*out_ << std::setw(10) << "Elapsed" << "\n";

		const Index width = 8 + 16 * dofNames.size() + 10;
		*out_ << tag() << std::string(std::max<Index>(30, width), '-') << "\n";

	}

} // namespace residuum::utils::logging::linear
