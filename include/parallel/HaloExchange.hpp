#ifndef RESIDUUM_PARALLEL_HALOEXCHANGE_HPP
#define RESIDUUM_PARALLEL_HALOEXCHANGE_HPP

#include <mpi.h>

#include "mesh/partition/PartitionData.hpp"

namespace residuum {
	namespace parallel {

		template<typename VectorT>
		class HaloExchange {
		public:
			HaloExchange(const PartitionData& partData, MPI_Comm comm);
			
			// reduce to host
			void reverseExchange(VectorT& vec) const;

			// broadcast ownder to ghosts
			void forwardExchange(VectorT& vec) const;

		private:
			struct SendRecvBuffer {};
			std::vector<SendRecvBuffer> sendBufs, recvBufs;
			MPI_Comm comm_;
		
		}; // class HaloExchange

	} // namespace parallel
} // namespace residuum

#endif
